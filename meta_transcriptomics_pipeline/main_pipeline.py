import argparse
import imp
import subprocess
import os
import operator
import time
import re
from get_lineage_info import get_lineage_info
#from main_pipeline.helpers import check_command_exists, check_fail, generate_temp_file
from helpers import check_command_exists, check_fail, generate_temp_file
from map_reads_to_contigs import map_reads_to_contigs
from merge_sams import merge_sams
from merge_contigs import merge_contigs
from remove_contaminants_control import remove_contaminants_control
from join_taxid_contigs import join_taxid_contigs
from map_reads_to_contigs import map_reads_to_contigs

def run_snap_single(index, in_path, out_path, threads):
    snap_contig_command = "snap-aligner" + " single " + index + " " + in_path +\
                    " -o " + out_path + " -t " + str(threads)
    new_command = subprocess.run(snap_contig_command, shell=True)
    if check_fail("snap-aligner", new_command, []) is True: return False
    
    return True

def run_diamond(index, in_path, out_path, threads, outfmt):
    diamond_command = "diamond" + " blastx --db " + index +\
                        " --query " + in_path + " --sensitive --max-target-seqs 1 --outfmt " + str(outfmt) +\
                        " --threads " + str(threads) +\
                        " --out " + out_path
    new_command = subprocess.run(diamond_command, shell=True)
    if check_fail("diamond", new_command, []) is True: return False

    return True

def fetch_contig_read_taxids(infile):
    assignments = {}
    with open(infile, "r") as r:
        for line in r:
            curr = line.split()
            assignments[curr[0]] = curr[2]

    return assignments

def compile_tpm_input_info(assignments, infile, outfile):
    if os.path.isfile(outfile):
        os.remove(outfile)
    wf = open(outfile, "w")
    with open(infile, "r") as r:
        for line in r:
            curr = line.split()
            if (curr[0] in assignments.keys()):
                name = curr[0]
                length = curr[1]
                count = curr[2]
                taxid = assignments[curr[0]]
                wf.write(name + "\t" + str(length) + "\t" + str(count) + "\t" + str(taxid) + "\n")

    wf.close()
    

def fetch_contig_taxids(infile):
    assignments = {}
    with open(infile, "r") as r:
        for line in r:
            curr = line.split()
            if (len(re.findall("^k[0-9]*", curr[0])) > 0):
                assignments[curr[0]] = curr[2]
           
    print("LEN " + str(len(assignments))) 
    return assignments

def assign_taxids_to_reads(assignments, infile, outfile):
    if os.path.isfile(outfile):
        os.remove(outfile)
    wf = open(outfile, "w")
    with open(infile, "r") as r:
        for line in r:
            #print(line)
            curr = line.split()
            name = curr[0]
            contig = curr[1]
            if (contig in assignments.keys()):
                wf.write(name + "\t" + str(assignments[contig]) + "\n")

    wf.close()

def assign_taxids_to_contigs_reads(assignments, infile, outfile):
    if os.path.isfile(outfile):
        os.remove(outfile)
    wf = open(outfile, "w")
    with open(infile, "r") as r:
        for line in r:
            #print(line)
            curr = line.split()
            name = curr[0]
            contig = curr[1]
            if (contig in assignments.keys()):
                wf.write(name + "\t" + str(assignments[contig]) + "\n")

    wf.close()              

# returns the number of reads assigned to each contig
def getContigReadCount(infile):
    counts = {}
    with open(infile, "r") as r:
        for line in r:
            curr = line.split()
            contig = curr[1]

            if contig in counts.keys():
                counts[contig] += 1
            else:
                counts[contig] = 1

    return counts

def countReads(infile, total_reads, outfile, contaminants):
    results = {}
    num_reads = 0
    print("total_reads: " + str(total_reads))
    with open(infile, "r") as f:
        for line in f:
            curr = line.split()
            taxid = curr[1]
            if taxid in results.keys():
                results[taxid] += 1
            else:
                results[taxid] = 1

            num_reads += 1

    for key in results.keys():
        results[key] = (results[key]/total_reads) * 100

    print(str(num_reads))
 
    unassigned_reads = total_reads - num_reads
    print("unassigned_reads: " + str(unassigned_reads))

    if (unassigned_reads != 0):
        results["Unknown"] = (unassigned_reads/total_reads) * 100
    
    sorted_results = dict( sorted(results.items(), key=operator.itemgetter(1),reverse=True))

    if os.path.isfile(outfile):
        os.remove(outfile)
    wf = open(outfile, "w")

    print(sorted_results)

    for key in sorted_results.keys():
        if (contaminants is not None and key not in contaminants):
            wf.write(key + "\t" + str(sorted_results[key]) + "\n")
        elif (contaminants is None):
            wf.write(key + "\t" + str(sorted_results[key]) + "\n")

    wf.close()

# infile must be in fastq
def getReadsLength(infile, outfile, contig_counts = None, readsTrue = False):

    wf = open(outfile, "a")
    with open(infile, "r") as f:
        divisible_int = 1
        curr_name = ""
        curr_seq = ""
        curr_length = -1
        for line in f:
            if (divisible_int%4 == 1):
                curr = line.split()
                curr_name = curr[0][1:].rstrip("\n")
            elif (divisible_int%4 == 2):
                curr_seq = line.rstrip("\n")
                curr_length = len(curr_seq)
            elif (divisible_int%4 == 3):
                # for reads
                if (readsTrue):
                    wf.write(curr_name + "\t" + str(curr_length) + "\t" + "1" + "\n")
                # for contigs
                else:
                    # assume contig_counts will never be NULL
                    # count = contig_counts[curr_name]
                    wf.write(curr_name + "\t" + str(curr_length) + "\t" + str(contig_counts[curr_name]) + "\n")
            else:
                divisible_int = 1
                continue

            divisible_int += 1

    wf.close()

def get_abundance(joined, total_reads, final_file, contaminants):
    final_res = {}
    taxids = {}
    lengths = {}
    counts = {}

    denom = 0.0

    with open(joined, "r") as f:
        for line in f:
            curr = line.split()
            name = curr[0]
            length = int(curr[1])
            count = int(curr[2])/total_reads
            taxid = int(curr[3])

            taxids[name] = taxid
            lengths[name] = length
            counts[name] = count
            denom += count
            
    # solution below adapted from RSEM
    # https://github.com/deweylab/RSEM/blob/e4dda70e90fb5eb9b831306f1c381f8bbf71ef0e/WriteResults.h
    # function adapted: calcExpressionValues
    # only difference is I did not use expected length, using the real length instead

    for key in lengths:
        counts[key] = counts[key]/denom

    fpkm = {}
    # calculate FPKM
    for key in lengths:
        fpkm[key] = (counts[key] * 1e9)/lengths[key]

    denom = 0.0
    for key in lengths:
        denom += fpkm[key]

    tpm = {}
    for key in lengths:
        tpm[key] = fpkm[key]/(denom * 1e6)

    for key in tpm:
        print(str(key) + "\t" + str(tpm[key]) + "\n")

    final_tpm = {}
    # calculating tpm but for the taxids
    for key in lengths:
        if taxids[key] in final_tpm.keys():
            final_tpm[taxids[key]] += tpm[key]
        else:
            final_tpm[taxids[key]] = tpm[key]

    if os.path.isfile(final_file):
        os.remove(final_file)
    wf = open(final_file, "w")
    final_tpm_sorted = dict( sorted(final_tpm.items(), key=operator.itemgetter(1),reverse=True))
    for key in final_tpm_sorted:
        if (key not in contaminants):
            wf.write(str(key) + "\t" + str(final_tpm_sorted[key]) + "\n")
    wf.close()

def run_pipeline(args: argparse.Namespace):

    if (len(args.input) > 2):
        print("Invalid file input")
        exit()

    input_files = []
    is_pe = False
    if (len(args.input) == 2):
        is_pe = True

    dirpath = args.dirpath

    ##################### FASTP ########################

    # generating temp files, to be removed later
    generated_files = []

    fastp_path = "fastp"
    fastp_out = []
    fastp_command = fastp_path

    # need to consider adapters, should we give the user a chance to add adatpers?

    if (is_pe is True):
        out1 = dirpath + "/qc_1.fastq"
        out2 = dirpath + "/qc_2.fastq"
        generated_files.append(out1)
        generated_files.append(out2)
        fastp_command = fastp_command +\
                        " --in1 " + args.input[0] +\
                        " --in2 " + args.input[1] +\
                        " --out1 " + out1 +\
                        " --out2 " + out2

        fastp_out.append(out1)
        fastp_out.append(out2)

    else:
        out = dirpath + "/qc.fastq"
        generated_files.append(out)
        fastp_command = fastp_command +\
                        " -i " + args.input[1] +\
                        " -o " + out

        fastp_out.append(out)

    
    fastp_command = fastp_command +\
                        " --qualified_quality_phred  " + args.qualified_quality_phred +\
                        " --unqualified_percent_limit " + args.unqualified_percent_limit +\
                        " --length_required " + args.length_required +\
                        " --detect_adapter_for_pe" +\
                        " --thread " + str(args.threads)

    start = time.time()
    new_command = subprocess.run(fastp_command, shell=True)
    if check_fail(fastp_path, new_command, generated_files) is True: return None
    end = time.time()
    print("fastp took: " + str(end - start))

    #################################### SORTMERNA ############################
    aligned = dirpath + "/aligned"
    other = dirpath + "/other"
    
    sortmerna_path = 'sortmerna'
    sortmerna_aligned = []
    sortmerna_other = []
    # best = 1 = all high candidate reference sequences will be searched for alignments
    sortmerna_command = sortmerna_path +\
                    " --ref " + args.sortmerna_index +\
                    " --aligned " + aligned +\
                    " --other " + other +\
                    " --fastx " +\
                    " --best 1 " +\
                    " --threads " + str(args.threads)

    if (is_pe == True):
        sortmerna_command = sortmerna_command +\
                    " --reads " + fastp_out[0] + " --reads " + fastp_out[1] +\
		            " --out2 TRUE " +\
		            " --paired_out TRUE"

        generated_files.append(other + "_fwd.fastq")
        generated_files.append(other + "_rev.fastq")
        sortmerna_aligned.append(aligned + "_fwd.fastq")
        sortmerna_aligned.append(aligned + "_rev.fastq")
        sortmerna_other.append(other + "_fwd.fastq")
        sortmerna_other.append(other + "_rev.fastq")

    else:
        sortmerna_command = sortmerna_command + " --reads " + fastp_out[0]

        generated_files.append(other + ".fastq")
        sortmerna_aligned.append(aligned + ".fastq")
        sortmerna_other.append(other + ".fastq")


    start = time.time()
    new_command = subprocess.run(sortmerna_command, shell=True)
    if check_fail(sortmerna_path, new_command, generated_files) is True: return None
    end = time.time()
    print("sortmerna took: " + str(end - start))
    #os.remove(aligned + "_fwd.fastq", aligned + ".log", aligned + "_rev.fastq")

    ############################# SNAP HUMAN #############################
    human_out = dirpath + "/snap_human_out.sam"
    snap_path = 'snap-aligner'
    snap_human_command = snap_path

    if (is_pe is True):
        snap_human_command = snap_human_command + " paired " + args.snap_human_index + " " + sortmerna_other[0] + " " + sortmerna_other[1]
    else:
        snap_human_command = snap_human_command + " single " + args.snap_human_index + " " + sortmerna_other[0]

    snap_human_command = snap_human_command +\
                        " -o " + human_out +\
                        " -t " + str(args.threads)

    start = time.time()
    new_command = subprocess.run(snap_human_command, shell=True)
    if check_fail(snap_path, new_command, [generated_files]) is True: return None
    end = time.time()
    print("human subtraction via snap took: " + str(end - start))

    human_subtract = []

    samtools_path = "samtools"
    samtools_human_subtract_command = samtools_path + " fastq  -f 4 -@ " + str(args.threads)

    if (is_pe is True):
        human_subtract.append(dirpath + "/human_subtract1.fastq")
        human_subtract.append(dirpath + "/human_subtract2.fastq")
        human_subtract.append(dirpath + "/human_spare.fastq")

        # retrieving only human reads
        samtools_human_subtract_command = samtools_human_subtract_command +\
                            " -1 " + human_subtract[0] +\
                            " -2 " + human_subtract[1] +\
                            " -s " + human_subtract[2] +\
                            " " + human_out

    else:
        human_subtract.append(dirpath + "/human_subtract.fastq")

        # retrieving only human reads
        samtools_human_subtract_command = samtools_human_subtract_command +\
                            " " + human_out +\
                            " > " + human_subtract[0]

    generated_files = generated_files + human_subtract

    new_command = subprocess.run(samtools_human_subtract_command, shell=True)
    if check_fail(samtools_path, new_command, []) is True: return None

    #################### MEGAHIT ###########################
    megahit_path = "megahit"
    contig_path = dirpath + "/megahit_out"
    megahit_command = megahit_path

    if (is_pe is True):
        megahit_command = megahit_command +\
                        " -1 " + human_subtract[0] +\
                        " -2 " + human_subtract[1]
    else:
        megahit_command = megahit_command + " -r " + human_subtract[0]

    megahit_command = megahit_command + " -o " + contig_path + " -t " + str(args.threads) # is an output directory 

    contigs = contig_path + "/final_contigs.fa"
    start = time.time()
    new_command = subprocess.run(megahit_command, shell=True)
    if check_fail(megahit_path, new_command, []) is True: return None
    end = time.time()
    print("assembly via megahit took: " + str(end - start))
    new_command = subprocess.run("mv " + contig_path + "/final.contigs.fa " + contigs, shell=True)

    contig_path = dirpath + "/megahit_out"
    contigs = contig_path + "/final_contigs.fa"
    human_subtract_1 = dirpath + "/human_subtract1.fastq"
    human_subtract_2 = dirpath + "/human_subtract2.fastq"
    samtools_path = "samtools"

    # need to convert file above from fa to fq, simply done using seqtk
    seqtk_path = "seqtk"
    new_contigs = contig_path + "/final_contigs.fq"
    seqtk_command = seqtk_path + " seq -F '#' " + contigs + " > " + new_contigs
    new_command = subprocess.run(seqtk_command, shell=True)
    if check_fail(seqtk_path, new_command, []) is True: return None

    # we can now align the contigs to the databases
    # lets do them sequentially for now
    snap_contigs = dirpath + "/snap_contigs_out"
    if run_snap_single(args.snap_index, new_contigs, snap_contigs, args.threads) == False: return None

    new_contigs_fa = contig_path + "/final_contigs.fa"
    blast_command = "blastn -query " + new_contigs_fa + " -db nt -out " + snap_contigs + " -outfmt 6"
    start = time.time()
    new_command = subprocess.run(blast_command, shell=True)
    if check_fail("blastn", new_command, []) is True: return None
    end = time.time()
    print("contig alignment against nt took: " + str(end - start))

    # now lets align reads
    diamond_contigs = dirpath + "/diamond_contigs_out"
    start = time.time()
    if run_diamond(args.diamond_index, new_contigs, diamond_contigs, args.threads, 6) == False: return None
    end = time.time()
    print("contig alignment against nr took: " + str(end - start))

    snap_contig_out, diamond_contig_out = merge_contigs(snap_contigs, diamond_contigs, dirpath)

    snap_contig_out = dirpath + "/nucl_alignments_contigs.txt"
    diamond_contig_out = dirpath + "/prot_alignments_contigs.txt"

    # we must retrieve the unaligned reads
    bbwrap_path = "bbwrap.sh"
    reads_mapped_to_contigs_file = dirpath + "/reads_mapped_to_contigs.sam"
    align_reads_to_contigs_cmd = bbwrap_path + " ref=" + contigs

    if (is_pe is True):
        align_reads_to_contigs_cmd = align_reads_to_contigs_cmd +\
                                    " in=" + human_subtract[0] +\
                                    " in2=" + human_subtract[1]
    else:
         align_reads_to_contigs_cmd = align_reads_to_contigs_cmd +\
                                    " in=" + human_subtract[0]

    align_reads_to_contigs_cmd = align_reads_to_contigs_cmd + " -out=" + reads_mapped_to_contigs_file  
    new_command = subprocess.run(align_reads_to_contigs_cmd, shell=True)
    if check_fail(bbwrap_path, new_command, []) is True: return None 

    # now lets retrieve the reads that did not align
    unaligned = []
    align_command = "samtools fastq -f 4 -@ " + str(args.threads)
    if (is_pe is True):
        unaligned.append(dirpath + "/new_fwd.fq")
        unaligned.append(dirpath + "/new_rev.fq")

        align_command = align_command + " -1 " + unaligned[0] +\
                        " -2 " + unaligned[1] + " " + reads_mapped_to_contigs_file

    else:
        unaligned.append(dirpath + "/new.fq")
        align_command = align_command + " " + reads_mapped_to_contigs_file + " > " +  unaligned[0]

    new_command = subprocess.run(align_command, shell=True)
    if check_fail(samtools_path, new_command, []) is True: return None 

    # running exact set of commands above but for unaligned reads this time

    snap_reads = dirpath + "/snap_reads_out.sam"
    snap_read_command = "snap-alinger"

    if (is_pe is True):
        snap_read_command = snap_read_command + " paired " + args.snap_index + " " + unaligned[0] + " " + unaligned[1]
    else: 
        snap_read_command = snap_read_command + " single " + args.snap_index + " " + unaligned[0]

    snap_read_command = snap_read_command + " -o " + snap_reads + " -I " + " -t " + str(args.threads)
    start = time.time()
    new_command = subprocess.run(snap_read_command, shell=True)
    if check_fail("snap-aligner", new_command, []) is True: return False
    end = time.time()
    print("read alignment against nt took: " + str(end - start))

    # now lets align reads
    # need to merge paired end reads first though

    diamond_input = ""
    if (is_pe is True):
        diamond_input = dirpath + "/merged_reads.fq"
        merge_command = "seqtk mergepe " + unaligned[0] + " " + unaligned[1] + " > " + diamond_input
        new_command = subprocess.run(merge_command, shell=True)
        if check_fail("seqtk mergepe", new_command, []) is True: return False 
    else: 
        diamond_input = unaligned[0]

    diamond_reads = dirpath + "/diamond_reads_out.sam"
    start = time.time()
    if run_diamond(args.diamond_index, diamond_input, diamond_reads, args.threads, 101) == False: return None
    end = time.time()
    print("read alignment against nr took: " + str(end - start))

    snap_reads_out, diamond_reads_out = merge_sams(snap_reads, diamond_reads, dirpath)

    #exit()

    snap_reads_out = dirpath + "/nucl_alignments_reads.txt"
    diamond_reads_out = dirpath + "/prot_alignments_reads.txt"
    
    # create 1 big file, for the reads and contigs mapped to their taxid, based on their accession_num
    contigs_reads_taxids_temp = dirpath + "/nucl_prot_taxids_temp.txt"
    if os.path.isfile(contigs_reads_taxids_temp):
        os.remove(contigs_reads_taxids_temp)
    if join_taxid_contigs(snap_contig_out, snap_reads_out, diamond_contig_out, diamond_reads_out, args.nucl_accession_taxid_mapping, args.prot_accession_taxid_mapping, contigs_reads_taxids_temp, dirpath) is False: return None
    contigs_reads_taxids_unsorted = dirpath + "/nucl_prot_taxids_unsorted.txt"
    subprocess.run("sed 's/ /\t/g' " + contigs_reads_taxids_temp + " > " + contigs_reads_taxids_unsorted, shell=True) # change space to tabs


    #exit()

    print("Checkpoint 1")

    # firstly lets count the reads
    num_reads_bytes = subprocess.run(['grep', '-c', '.*', human_subtract_1], capture_output=True)
    num_reads_str = num_reads_bytes.stdout.decode('utf-8')
    num_reads = int(num_reads_str.replace('\n', ''))/4 # finally in int format, dividing by 4 because its in fastq format

    # map reads to their contigs
    mapped_reads_unsorted = dirpath + "/reads_mapped_to_contigs_unsorted.txt"
    map_reads_to_contigs(reads_mapped_to_contigs_file, mapped_reads_unsorted, dirpath)

    # #exit()

    # now lets map the contig taxids to their reads
    # firstly retrieve all the dna/protein aligned reads

    reads_taxids_temp = dirpath + "/all_reads_taxids_temp.txt"
    contigs_reads_taxids = dirpath + "/nucl_prot_taxids.txt"
    subprocess.run("LC_COLLATE=C sort -k 1 " + contigs_reads_taxids_unsorted + " > " + contigs_reads_taxids, shell=True)
    
    print("Checkpoint 2")

    mapped_reads = dirpath + "/reads_mapped_to_contigs.txt"
    subprocess.run("LC_COLLATE=C sort -k 2 " + mapped_reads_unsorted + " > " + mapped_reads, shell=True)
    assignments = fetch_contig_taxids(contigs_reads_taxids)
    assign_taxids_to_reads(assignments, mapped_reads, reads_taxids_temp)
    reads_taxids = dirpath + "/all_reads_taxids.txt"
    subprocess.run("sed 's/ /\t/g' " + reads_taxids_temp + " > " + reads_taxids, shell=True) # change space to tabs

    # we need to now combine the contig aligned reads to the non contig aligned reads
    subprocess.run("awk '$1 !~ /^k[0-9]*/ {print $1\"\t\"$3}' " + contigs_reads_taxids + " >> " + reads_taxids, shell=True)


    # before proceeeding any further, lets remove the contaminants, simply remove the taxids
    contaminant_removal = True

    print("Hey?")

    if args.kraken_db is None:
        print("Path to kraken index is not provided, skipping contamination removal")
        contaminant_removal = False

    if args.control_sequences is None or args.other_sequences is None or len(args.control_sequences + args.other_sequences) == 0 and contaminant_removal == True:
        print("No control/additional sample files have been provided, skipping contamination removal")
        contaminant_removal = False

    print("Starting decontamination")

    contaminants = []
    rcf_out = dirpath + "/rcf_out.txt"
    if contaminant_removal is True:
        others = args.other_sequences
        others.append(args.inp1)
        others.append(args.inp2)
        contaminants = remove_contaminants_control(args.control_sequences, others, args.kraken_db, rcf_out, args.taxdump_location, args.threads, dirpath)


    print("Finished decontam")

    # now we can calculate read count method
    readCountsOutfile = dirpath + "/readCountsOut.txt"
    countReads(reads_taxids, num_reads, readCountsOutfile, contaminants)

    # now lets do abundance calculations via the tpm method
    # firstly get the length of the contigs
    # then get the number of reads aligned to them

    # firstly lets get the number of reads assigned to each contig
    contig_counts = getContigReadCount(mapped_reads)

    contig_unaligned_read_counts_temp = dirpath + "/contig_unaligned_read_counts_temp.txt"
    if os.path.isfile(contig_unaligned_read_counts_temp):
        os.remove(contig_unaligned_read_counts_temp)
    getReadsLength(new_contigs, contig_unaligned_read_counts_temp, contig_counts, False)

    # then repeat process of unaligned reads
    getReadsLength(new_fwd, contig_unaligned_read_counts_temp, None, True)

    # lets join contig_unaligned_read_counts with their taxid
    contig_unaligned_read_counts_temp2 = dirpath + "/contig_unaligned_read_counts_temp2.txt"
    assignments = fetch_contig_read_taxids(contigs_reads_taxids)
     

    contig_unaligned_read_counts = dirpath + "/contig_unaligned_read_counts.txt"
    compile_tpm_input_info(assignments, contig_unaligned_read_counts_temp, contig_unaligned_read_counts_temp2)
    subprocess.run("sed 's/ /\t/g' " + contig_unaligned_read_counts_temp2 + " > " + contig_unaligned_read_counts, shell=True) # change space to tabs

    # now we can finally calculate TPM/FPKM
    tpm_abundance_file = dirpath + "/tpm_fpkm.txt"
    get_abundance(contig_unaligned_read_counts, num_reads, tpm_abundance_file, contaminants)

    #exit()

    print("DONE!!!")

    print("In the white room")

    # now lets get plot the abundances as krona charts
    # need to find lineages first
    readAbundances = dirpath + "/readAbundances.txt"
    get_lineage_info(readCountsOutfile, readAbundances, args.taxdump_location)
    readAbundancesKrona = dirpath + "/readAbundancesKrona.html"
    subprocess.run("ktImportText " + readAbundances + " -o " + readAbundancesKrona, shell=True)

    tpmAbundances = dirpath + "/tpmAbundances.txt"
    get_lineage_info(tpm_abundance_file, tpmAbundances, args.taxdump_location)
    tpmAbundancesKrona = dirpath + "/tpmAbundancesKrona.html"
    subprocess.run("ktImportText " + tpmAbundances + " -o " + tpmAbundancesKrona, shell=True)

    print("DONE!!!")
