import argparse
import subprocess
import time
import os
import re 
import operator 
from get_lineage_info import get_lineage_info
#from main_pipeline.helpers import check_command_exists, check_fail, generate_temp_file
from map_reads_to_contigs import map_reads_to_contigs
from merge_blast_outputs import merge_blast_outputs
from remove_contaminants_control import remove_contaminants_control
from join_taxid_contigs import join_taxid_contigs
from map_reads_to_contigs import map_reads_to_contigs
from match_scores import match_scores
from filter_results import filter_result, get_filtered_taxids
from get_abundance import get_abundance

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
                    if curr_name in contig_counts:
                        wf.write(curr_name + "\t" + str(curr_length) + "\t" + str(contig_counts[curr_name]) + "\n")
            else:
                divisible_int = 1
                continue

            divisible_int += 1

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

def finalisation(args: argparse.Namespace):
    dirpath = args.dirpath
    nt_out = dirpath + "/nucl_alignments_contigs.txt"
    nr_out = dirpath + "/prot_alignments_contigs.txt"
    nt_dummy = dirpath + "/nt_dummy.txt"
    nr_dummy = dirpath + "/nr_dummy.txt"
    new_fwd = dirpath + "/new_fwd.fq"
    human_subtract_1 = dirpath + "/human_subtract1.fastq"
    reads_mapped_to_contigs_file = dirpath + "/reads_mapped_to_contigs.sam"
    new_contigs = contig_path + "/final_contigs.fq"
    contig_path = dirpath + "/megahit_out"
    num_reads_bytes = subprocess.run(['grep', '-c', '.*', human_subtract_1], capture_output=True)
    num_reads_str = num_reads_bytes.stdout.decode('utf-8')
    num_reads = int(num_reads_str.replace('\n', ''))/4 # finally in int format, dividing by 4 because its in fastq format

    subprocess.run("touch " + nt_dummy, shell=True)
    subprocess.run("touch " + nr_dummy, shell=True)
    
    # create 1 big file, for the reads and contigs mapped to their taxid, based on their accession_num
    start = time.time()
    contigs_reads_taxids_temp = dirpath + "/nucl_prot_taxids_temp.txt"
    if os.path.isfile(contigs_reads_taxids_temp):
        os.remove(contigs_reads_taxids_temp)
    if join_taxid_contigs(nt_out, nt_dummy, nr_out, nr_dummy, args.nucl_accession_taxid_mapping, args.prot_accession_taxid_mapping, contigs_reads_taxids_temp, dirpath) is False: return None
    contigs_reads_taxids_unsorted = dirpath + "/nucl_prot_taxids_unsorted.txt"
    subprocess.run("sed 's/ /\t/g' " + contigs_reads_taxids_temp + " > " + contigs_reads_taxids_unsorted, shell=True) # change space to tabs
    end = time.time()
    print("taxid identification via accessions took: " + str(end - start))

    # map reads to their contigs
    mapped_reads_unsorted = dirpath + "/reads_mapped_to_contigs_unsorted.txt"
    map_reads_to_contigs(reads_mapped_to_contigs_file, mapped_reads_unsorted, dirpath)

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
        start = time.time()
        contaminants = remove_contaminants_control(args.control_sequences, args.other_sequences, args.kraken_db, rcf_out, args.taxdump_location, args.threads, dirpath)
        end = time.time()
        print("contaminant identification (kraken + recentrifuge) took: " + str(end - start))

    print("Finished decontamination")

    # now we can calculate read count method
    readCountsOutfile = dirpath + "/readCountsOut.txt"
    countReads(reads_taxids, num_reads, readCountsOutfile, contaminants)

    to_filter_out = filter_result(dirpath + "/temp_out_2", args.pid_filter, args.evalue_filter, args.bitscore_filter)
    readCountsFiltered = dirpath + "/readCountsFiltered.txt"
    get_filtered_taxids(to_filter_out, readCountsOutfile, readCountsFiltered)

    # now lets get plot the abundances as krona charts
    # need to find lineages first
    readAbundances = dirpath + "/readAbundances.txt"
    get_lineage_info(readCountsOutfile, readAbundances, args.taxdump_location)
    readAbundancesKrona = dirpath + "/readAbundancesKrona.html"
    subprocess.run("ktImportText " + readAbundances + " -o " + readAbundancesKrona, shell=True)
    subprocess.run("ImportText.pl " + readAbundances + " -o " + readAbundancesKrona + " -fil " + taxid_scores, shell=True)

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

    log_path = dirpath + "/megahit_out/log"
    n50 = subprocess.check_output('tail -n 2 ' + log_path + ' | head -n 1', shell=True).decode('utf-8').strip('\n').split(' ')[-2:-1][0]

    # now we can finally calculate TPM/FPKM
    tpm_abundance_file = dirpath + "/tpm_fpkm.txt"
    get_abundance(contig_unaligned_read_counts, num_reads, n50, tpm_abundance_file, contaminants)

    tpmAbundances = dirpath + "/tpmAbundances.txt"
    get_lineage_info(tpm_abundance_file, tpmAbundances, args.taxdump_location)
    tpmAbundancesKrona = dirpath + "/tpmAbundancesKrona.html"
    subprocess.run("ImportText.pl " + tpmAbundances + " -o " + tpmAbundancesKrona  + " -fil " + taxid_scores, shell=True)