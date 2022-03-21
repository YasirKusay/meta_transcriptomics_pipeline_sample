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
from paf2blast6 import paf2blast6
from filter_files import filter_files
from match_scores import match_scores
from filter_results import filter_result, get_filtered_taxids


def run_snap_single(index, in_path, out_path, threads):
    snap_contig_command = "snap-aligner" + " single " + index + " " + in_path +\
                    " -o " + out_path + " -t " + str(threads)
    new_command = subprocess.run(snap_contig_command, shell=True)
    if check_fail("snap-aligner", new_command, []) is True: return False
    
    return True

def run_diamond(index, in_path, out_path, threads, outfmt):
    diamond_command = "diamond" + " blastx --db " + index +\
                        " --query " + in_path + " --sensitive --max-target-seqs 1 --outfmt " + str(outfmt) +\
                        " --masking 0 -c 1 -s 1 -b 0.75" +\
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

def process_fast_mode_output(kraken_out, outfile, total_reads):
    results = {}
    num_reads = 0
    with open(kraken_out, "r") as f:   
        for line in f:
            curr = line.split()
            taxid = curr[2]
            if taxid == 0:
                taxid == "Unknown"
            if taxid in results.keys():
                results[taxid] += 1
            else:
                results[taxid] = 1

            num_reads += 1

    for key in results.keys():
        results[key] = (results[key]/total_reads) * 100

    sorted_results = dict( sorted(results.items(), key=operator.itemgetter(1),reverse=True))
    wf = open(outfile, "w")

    print(sorted_results)

    for key in sorted_results.keys():
        wf.write(key + "\t" + str(sorted_results[key]) + "\n")
        wf.write(key + "\t" + str(sorted_results[key]) + "\n")

    wf.close()

def run_pipeline(args: argparse.Namespace):

    # setting up a temporary dir for all intermediate files generated
    #dirpath = tempfile.mkdtemp()

    # lets instead check if pipeline works ok :)
    dirpath = args.dirpath

    ##################### FASTP ########################

    # generating temp files, to be removed later
    #out1 = "/srv/scratch/z5215055/HONS/res/out1.fastq"
    #out2 = "/srv/scratch/z5215055/HONS/res/out2.fastq"
    generated_files = []

    out1 = dirpath + "/qc_1.fastq"
    out2 = dirpath + "/qc_2.fastq"
    generated_files.append(out1)
    generated_files.append(out2)

    fastp_path = "fastp"
    fastp_command = fastp_path +\
                    " --in1 " + args.inp1 +\
                    " --in2 " + args.inp2 +\
                    " --out1 " + out1 +\
                    " --out2 " + out2 +\
                    " --qualified_quality_phred  " + args.qualified_quality_phred +\
                    " --unqualified_percent_limit " + args.unqualified_percent_limit +\
                    " --length_required " + args.length_required +\
                    " --detect_adapter_for_pe" +\
                    " --thread " + str(args.threads)
                    # need to consider adapters, should we give the user a chance to add adatpers?

    start = time.time()
    #new_command = subprocess.run(fastp_command, shell=True)
    #if check_fail(fastp_path, new_command, [out1, out2]) is True: return None
    end = time.time()
    print("fastp took: " + str(end - start))

    #################################### SORTMERNA ############################
    aligned = dirpath + "/aligned"
    other = dirpath + "/other"
    
    sortmerna_path = 'sortmerna'
    sortmerna_command = sortmerna_path +\
                    " --ref " + args.sortmerna_index +\
                    " --aligned " + aligned +\
                    " --other " + other +\
                    " --fastx " +\
                    " --reads " + out1 + " --reads " + out2 +\
                    " --threads " + str(args.threads) +\
		            " --out2 TRUE " +\
		            " --paired_out TRUE" +\
                    " --best 1 " # 1 = all high candidate reference sequences will be searched for alignments

    start = time.time()
    #new_command = subprocess.run(sortmerna_command, shell=True)
    #if check_fail(sortmerna_path, new_command, [out1, out2, aligned + "_fwd.fastq", aligned + "_rev.fastq",  other + "_fwd.fastq", other + "_rev.fastq"]) is True: return None
    end = time.time()
    print("sortmerna took: " + str(end - start))
    #os.remove(aligned + "_fwd.fastq", aligned + ".log", aligned + "_rev.fastq")
    generated_files.append(other + "_fwd.fastq")
    generated_files.append(other + "_rev.fastq")

    ############################# SNAP HUMAN #############################
    human_out = dirpath + "/snap_human_out.sam"
    snap_path = 'snap-aligner'
    snap_human_command = snap_path + " paired " + args.snap_human_index + " " + other + "_fwd.fastq " + other + "_rev.fastq " +\
                    " -o " + human_out + " -t " + str(args.threads) + " -I "
    start = time.time()
    #new_command = subprocess.run(snap_human_command, shell=True)
    #if check_fail(snap_path, new_command, [generated_files]) is True: return None
    end = time.time()
    print("human subtraction via snap took: " + str(end - start))

    human_subtract_1 = dirpath + "/human_subtract1.fastq"
    human_subtract_2 = dirpath + "/human_subtract2.fastq"
    human_spare = dirpath + "/human_spare.fastq"

    # retrieving only human reads
    samtools_path = "samtools"
    samtools_human_subtract_command = samtools_path + " fastq  -f 4 -@ " + str(args.threads) +\
                        " -1 " + human_subtract_1 +\
                        " -2 " + human_subtract_2 +\
                        " -s " + human_spare + " " +\
                        human_out

    #new_command = subprocess.run(samtools_human_subtract_command, shell=True)
    #if check_fail(samtools_path, new_command, []) is True: return None
    generated_files.append(human_subtract_1 + ".fastq")
    generated_files.append(human_subtract_2 + ".fastq")
    generated_files.append(human_spare + ".fastq")

    # QUICK ALIGNMENT, JUST ALIGN REMAINING READS USING KRAKEN AGAINST KRAKEN_PLUS
    num_reads_bytes = subprocess.run(['grep', '-c', '.*', human_subtract_1], capture_output=True)
    num_reads_str = num_reads_bytes.stdout.decode('utf-8')
    num_reads = int(num_reads_str.replace('\n', ''))/4 # finally in int format, dividing by 4 because its in fastq format
    fast_mode_output = dirpath + "/fast_mode_output"
    kraken_command = "kraken2 --db " + args.kraken_db + " --threads " + str(args.threads) +\
                        " --output " + fast_mode_output + " --paired " + human_subtract_1 + " " + human_subtract_2
    #new_command = subprocess.run(kraken_command, shell=True)
    #if check_fail("kraken", new_command, []) is True: return None
    kraken_res_out = dirpath + "/kraken_res_out"
    #process_fast_mode_output(fast_mode_output, kraken_res_out, num_reads)
    fastAbundances = dirpath + "/fastAbundances.txt"
    #get_lineage_info(kraken_res_out, fastAbundances, args.taxdump_location)
    #fastAbundancesKrona = dirpath + "/fastAbundancesKrona.html"
    #subprocess.run("ImportText.pl " + fastAbundances + " -o " + fastAbundancesKrona, shell=True)
    #exit()

    #################### MEGAHIT ###########################
    megahit_path = "megahit"
    contig_path = dirpath + "/megahit_out"
    megahit_command = megahit_path + " -1 " + human_subtract_1 +\
                        " -2 " + human_subtract_2 +\
                        " -o " + contig_path + " -t " + str(args.threads) # is an output directory 
    contigs = contig_path + "/final_contigs.fa"
    start = time.time()
    #new_command = subprocess.run(megahit_command, shell=True)
    #if check_fail(megahit_path, new_command, []) is True: return None
    end = time.time()
    print("assembly via megahit took: " + str(end - start))
    new_command = subprocess.run("mv " + contig_path + "/final.contigs.fa " + contigs, shell=True)

    # we must retrieve the unaligned reads
    bbwrap_path = "bbwrap.sh"
    reads_mapped_to_contigs_file = dirpath + "/reads_mapped_to_contigs.sam"
    align_reads_to_contigs_cmd = bbwrap_path + " ref=" + contigs +\
                                " in=" + human_subtract_1 +\
                                " in2=" + human_subtract_2 +\
                                " -out=" + reads_mapped_to_contigs_file  
    #new_command = subprocess.run(align_reads_to_contigs_cmd, shell=True)
    #if check_fail(bbwrap_path, new_command, []) is True: return None 

    # now lets retrieve the reads that did not align
    new_fwd = dirpath + "/new_fwd.fq"
    new_rev = dirpath + "/new_rev.fq"

    align_command = "samtools fastq -f4 -1 " + new_fwd +\
                    " -2 " + new_rev + " " + reads_mapped_to_contigs_file
    #new_command = subprocess.run(align_command, shell=True)
    #if check_fail(samtools_path, new_command, []) is True: return None 

    smaller_1 = dirpath + "/smaller_1.fq"
    smaller_2 = dirpath + "/smaller_2.fq"
    bigger_1 = dirpath + "/bigger_1.fq"
    bigger_2 = dirpath + "/bigger_2.fq"

    filter_files(new_fwd, bigger_1, smaller_1)
    filter_files(new_rev, bigger_2, smaller_2)

    #time minimap2 -a -H -t 32 /data/bio/ncbi/current/nt.mmi /srv/scratch/z5215055/HONS/samples/csf/megahit_out/final_contigs.fa > /srv/scratch/z5215055/HONS/minimap2/output.paf
    # minimap2 -ax sr ref.fa read1.fq read2.fq > aln.sam  

    # need to convert file above from fa to fq, simply done using seqtk
    seqtk_path = "seqtk"
    new_contigs = contig_path + "/final_contigs.fq"
    seqtk_command = seqtk_path + " seq -F '#' " + contigs + " > " + new_contigs
    #new_command = subprocess.run(seqtk_command, shell=True)
    #if check_fail(seqtk_path, new_command, []) is True: return None

    nt_combined_file = dirpath + "/nt_combined_file"

    # aligning contigs + longer reads, need to merge big reads first
    start = time.time()
    minimap2_path = "minimap2"
    minimap2_contig_out = dirpath + "/minimap2_contig_out.paf"
    minimap2_command = minimap2_path + " -c -t " + str(args.threads) + " --split-prefix " + dirpath + " " + args.minimap2_index + " " + new_contigs + " > " + minimap2_contig_out
    #new_command = subprocess.run(minimap2_command, shell=True)
    #if check_fail(minimap2_path, new_command, []) is True: return None


    end = time.time()
    print("minimap2 alignment 1 against nt took: " + str(end - start))

    start = time.time()
    paf2blast6(minimap2_contig_out, dirpath)
    paf2blast_out = dirpath + "/minimap2_contig_out_frompaf.m8"

    new_command = subprocess.run("cat " + paf2blast_out + " > " + nt_combined_file, shell=True)

    minimap2_long_reads_out = dirpath + "/minimap2_lr_out.paf"
    minimap2_command = minimap2_path + " -c -x sr -t " + str(args.threads) + " --split-prefix " + dirpath + " " + args.minimap2_index + " " + bigger_1 + " " + bigger_2 + " > " + minimap2_long_reads_out
    #new_command = subprocess.run(minimap2_command, shell=True)
    #if check_fail(minimap2_path, new_command, []) is True: return None

    end = time.time()
    print("minimap2 alignment 2 against nt took: " + str(end - start))

    paf2blast6(minimap2_long_reads_out, dirpath)
    paf2blast_out = dirpath + "/minimap2_lr_out_frompaf.m8"

    new_command = subprocess.run("cat " + paf2blast_out + " >> " + nt_combined_file, shell=True)

    # now lets blast the small reads, (merge them first)

    merged_short_pe = dirpath + "/merged_short_reads.fq"
    merge_command = "seqtk mergepe " + smaller_1 + " " + smaller_2 + " > " + merged_short_pe
    new_command = subprocess.run(merge_command, shell=True)
    if check_fail("seqtk mergepe", new_command, []) is True: return False 

    merged_short_pe_fa = dirpath + "/merged_short_reads.fa"
    seqtk_command = seqtk_path + " seq -a " + merged_short_pe + " > " + merged_short_pe_fa
    new_command = subprocess.run(seqtk_command, shell=True)
    if check_fail("seqtk", new_command, []) is True: return False 

    start = time.time()

    blast_path = "blastn -task megablast"
    blast_short_out = dirpath + "/blast_short_out"
    blast_command = blast_path + " -query " + merged_short_pe_fa + " -db nt " + " -out " + blast_short_out + " -outfmt 6 -max_target_seqs 20 -num_threads " + str(args.threads)
    #new_command = subprocess.run(blast_command, shell=True)
    #if check_fail(blast_path, new_command, []) is True: return False  

    end = time.time()
    print("blast alignment against nt took: " + str(end - start))


    new_command = subprocess.run("cat " + blast_short_out + " >> " + nt_combined_file, shell=True)

    # now repeat process for diamond

    # now lets align reads
    # need to merge paired end reads first though
    merged_pe = dirpath + "/merged_reads.fq"
    merge_command = "seqtk mergepe " + new_fwd + " " + new_rev + " > " + merged_pe
    new_command = subprocess.run(merge_command, shell=True)
    if check_fail("seqtk mergepe", new_command, []) is True: return False 

    diamonc_combined_file = dirpath + "/diamond_combined_file.fq"
    subprocess.run("cat " + new_contigs + " > " + diamonc_combined_file, shell=True)
    subprocess.run("cat " + merged_pe + " >> " + diamonc_combined_file, shell=True)

    nr_combined_file = dirpath + "/nr_combined_file"
    start = time.time()
    #if run_diamond(args.diamond_index, diamonc_combined_file, nr_combined_file, args.threads, 6) == False: return None
    end = time.time()
    print("contig alignment against nr took: " + str(end - start))

    nt_out, nr_out = merge_contigs(nt_combined_file, nr_combined_file, dirpath)

    nt_out = dirpath + "/nucl_alignments_contigs.txt"
    nr_out = dirpath + "/prot_alignments_contigs.txt"
    nt_dummy = dirpath + "/nt_dummy.txt"
    nr_dummy = dirpath + "/nr_dummy.txt"

    subprocess.run("touch " + nt_dummy, shell=True)
    subprocess.run("touch " + nr_dummy, shell=True)
    
    # create 1 big file, for the reads and contigs mapped to their taxid, based on their accession_num
    contigs_reads_taxids_temp = dirpath + "/nucl_prot_taxids_temp.txt"
    if os.path.isfile(contigs_reads_taxids_temp):
        os.remove(contigs_reads_taxids_temp)
    if join_taxid_contigs(nt_out, nt_dummy, nr_out, nr_dummy, args.nucl_accession_taxid_mapping, args.prot_accession_taxid_mapping, contigs_reads_taxids_temp, dirpath) is False: return None
    contigs_reads_taxids_unsorted = dirpath + "/nucl_prot_taxids_unsorted.txt"
    subprocess.run("sed 's/ /\t/g' " + contigs_reads_taxids_temp + " > " + contigs_reads_taxids_unsorted, shell=True) # change space to tabs


    #exit()

    print("Checkpoint 1")

    # map reads to their contigs
    mapped_reads_unsorted = dirpath + "/reads_mapped_to_contigs_unsorted.txt"
    map_reads_to_contigs(reads_mapped_to_contigs_file, mapped_reads_unsorted, dirpath)

    # #exit()

    # now lets map the contig taxids to their reads
    # firstly retrieve all the dna/protein aligned reads

    reads_taxids_temp = dirpath + "/all_reads_taxids_temp.txt"
    contigs_reads_taxids = dirpath + "/nucl_prot_taxids.txt"
    subprocess.run("LC_COLLATE=C sort -k 1 " + contigs_reads_taxids_unsorted + " > " + contigs_reads_taxids, shell=True)

    taxid_scores = dirpath + "/taxid_scores"
    combined_nt_nr_unsorted = dirpath + "/combined_nt_nr_unsorted"
    subprocess.run("cat " + nt_out + " > " + combined_nt_nr_unsorted, shell=True)
    subprocess.run("cat " + nr_out + " >> " + combined_nt_nr_unsorted, shell=True)
    combined_nt_nr = dirpath + "/combined_nt_nr"
    subprocess.run("LC_COLLATE=C sort -k 1 " + combined_nt_nr_unsorted + " > " + combined_nt_nr, shell=True)
    match_scores(contigs_reads_taxids, combined_nt_nr, dirpath, taxid_scores, args.taxdump_location)
    
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
        #others = args.other_sequences
        #others.append(args.inp1)
        #others.append(args.inp2)
        contaminants = remove_contaminants_control(args.control_sequences, args.other_sequences, args.kraken_db, rcf_out, args.taxdump_location, args.threads, dirpath)

    #exit()

    print("Finished decontam")

    # now we can calculate read count method
    readCountsOutfile = dirpath + "/readCountsOut.txt"
    countReads(reads_taxids, num_reads, readCountsOutfile, contaminants)

    to_filter_out = filter_result(dirpath + "/temp_out_2", args.pid_filter, args.evalue_filter, args.bitscore_filter)
    readCountsFiltered = dirpath + "/readCountsFiltered.txt"
    get_filtered_taxids(to_filter_out, readCountsOutfile, readCountsFiltered)

    readAbundances = dirpath + "/readAbundances.txt"
    get_lineage_info(readCountsFiltered, readAbundances, args.taxdump_location)
    readAbundancesKrona = dirpath + "/readAbundancesKrona.html"
    #subprocess.run("ktImportText " + readAbundances + " -o " + readAbundancesKrona, shell=True)
    subprocess.run("ImportText.pl " + readAbundances + " -o " + readAbundancesKrona + " -fil " + taxid_scores, shell=True)

    #exit()

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
    #readAbundances = dirpath + "/readAbundances.txt"
    #get_lineage_info(readCountsOutfile, readAbundances, args.taxdump_location)
    #readAbundancesKrona = dirpath + "/readAbundancesKrona.html"
    #subprocess.run("ktImportText " + readAbundances + " -o " + readAbundancesKrona, shell=True)

    tpmAbundances = dirpath + "/tpmAbundances.txt"
    get_lineage_info(tpm_abundance_file, tpmAbundances, args.taxdump_location)
    tpmAbundancesKrona = dirpath + "/tpmAbundancesKrona.html"
    subprocess.run("ImportText.pl " + tpmAbundances + " -o " + tpmAbundancesKrona  + " -fil " + taxid_scores, shell=True)

    print("DONE!!!")
