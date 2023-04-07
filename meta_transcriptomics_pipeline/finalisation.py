import argparse
import subprocess
import time
import os
import re 
import operator 
from meta_transcriptomics_pipeline.get_lineage_info import get_lineage_info
from meta_transcriptomics_pipeline.get_best_sam_hits import get_best_sam_hits
from meta_transcriptomics_pipeline.get_best_blast_hits import get_best_blast_hits
from meta_transcriptomics_pipeline.remove_contaminants_control import remove_contaminants_control
from meta_transcriptomics_pipeline.join_seq_to_taxid import join_seq_to_taxid
from meta_transcriptomics_pipeline.match_scores import match_scores
from meta_transcriptomics_pipeline.get_abundance import get_abundance
from meta_transcriptomics_pipeline.get_abundance import get_abundance
from meta_transcriptomics_pipeline.get_abundance import get_abundance
from meta_transcriptomics_pipeline.count_num_seqs import countNumLines

def fetch_taxids(infile):
    taxids = []
    with open(infile, "r") as f:
        for line in f:
            curr = line.split("\t")
            if curr[2].strip() not in taxids:
                taxids.append(curr[2].strip())

    return taxids

def fetch_contig_read_taxids(infile):
    assignments = {}
    with open(infile, "r") as r:
        for line in r:
            curr = line.split()
            assignments[curr[0]] = curr[2]

    return assignments

def compile_tpm_input_info(assignments, infile, outfile):
    wf = open(outfile, "w")
    with open(infile, "r") as r:
        for line in r:
            curr = line.split()
            if (curr[0] in assignments.keys()):
                name = curr[0]
                length = curr[1]
                count = curr[2].strip()
                taxid = assignments[curr[0]]
                wf.write(name + "\t" + str(length) + "\t" + str(count) + "\t" + str(taxid) + "\n")

    wf.close()
    

def fetch_contig_taxids(infile):
    assignments = {}
    with open(infile, "r") as r:
        for line in r:
            curr = line.split()
            if (len(re.findall("^k[0-9]*", curr[0])) > 0): # syntax for megahit contig names
                assignments[curr[0]] = curr[2].strip()

    return assignments

def assign_taxids_to_assembled_reads(assignments, infile, outfile):
    wf = open(outfile, "w")
    with open(infile, "r") as r:
        for line in r:
            #print(line)
            curr = line.split()
            name = curr[0]
            contig = curr[1]
            if (contig in assignments.keys()):
                wf.write(name + "\t" + str(assignments[contig].strip()) + "\n")

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

def countReads(infile, total_reads, outfile, contaminants, bad_taxids):
    results = {}
    num_reads = 0
    with open(infile, "r") as f:
        for line in f:
            curr = line.split()
            taxid = curr[1].strip()
            if taxid in results.keys():
                results[taxid] += 1
            else:
                results[taxid] = 1

            num_reads += 1

    for key in results.keys():
        results[key] = (results[key]/total_reads) * 100
 
    unassigned_reads = total_reads - num_reads
    print("unassigned_reads: " + str(unassigned_reads))

    # lets exclude unknown from final output
    # if (unassigned_reads != 0):
    #    results["Unknown"] = (unassigned_reads/total_reads) * 100
    
    sorted_results = dict( sorted(results.items(), key=operator.itemgetter(1),reverse=True))

    wf = open(outfile, "w")

    for taxid in sorted_results.keys():
        if (contaminants is not None and taxid not in contaminants and taxid not in bad_taxids):
            wf.write(taxid + "\t" + str(sorted_results[taxid]) + "\n")
        elif (contaminants is None and taxid not in bad_taxids):
            wf.write(taxid + "\t" + str(sorted_results[taxid]) + "\n")

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
                curr_name = curr[0][1:].rstrip("\n") # getting rid of statistics if its a contig
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

def finalisation(args: argparse.Namespace):
    dirpath = args.dirpath

    if dirpath[-1] == "/":
        dirpath = dirpath[0:-1]

    analysis_path = dirpath + "/analysis"

    if os.path.exists(analysis_path) is False:
        os.mkdir(analysis_path)

    final_plots_path = dirpath + "/final_plots"

    if os.path.exists(final_plots_path) is False:
        os.mkdir(final_plots_path)

    fullyQc1 = dirpath + "/preprocessing/fullyQc_fwd.fq"

    num_reads_bytes = subprocess.check_output(['grep', '-c', '.*', fullyQc1])
    num_reads_str = num_reads_bytes.decode('utf-8')
    num_reads = int(num_reads_str.replace('\n', ''))/4 # finally in int format, dividing by 4 because its in fastq format

    reads_mapped_to_contigs_file = dirpath + "/preprocessing/reads_mapped_to_contigs.sam"
    nt_alignments_file = dirpath + "/alignments/nt_alignments_file.tsv"
    nr_alignments_file = dirpath + "/alignments/nr_alignments_file.tsv"
    contigs_fq = dirpath + "/megahit_out/final_contigs.fq"

    # paf2blast_out_contigs = dirpath + "/minimap2_contig_out_frompaf.m8"
    # new_command = subprocess.run("cat " + paf2blast_out_contigs + " >> " + nt_alignments_file, shell=True)
    # paf2blast_out_lr = dirpath + "/minimap2_lr_out_frompaf.m8"
    # new_command = subprocess.run("cat " + paf2blast_out_lr + " >> " + nr_alignments_file, shell=True)

    best_nt_scores = analysis_path + "/best_nt_scores.tsv"
    best_nr_scores = analysis_path + "/best_nr_scores.tsv"
    get_best_blast_hits(nt_alignments_file, nr_alignments_file, analysis_path, best_nt_scores, best_nr_scores)

    nucl_accession_taxid_mapping_files = []
    prot_accession_taxid_mapping_files = []

    for currFile in os.listdir(args.nucl_prot_accession_taxid_mapping_files_loc):
        toAppend = ""
        if args.nucl_prot_accession_taxid_mapping_files_loc[-1] != "/":
            toAppend = args.nucl_prot_accession_taxid_mapping_files_loc + "/" + currFile
        else:
            toAppend = args.nucl_prot_accession_taxid_mapping_files_loc + currFile
        if os.path.isfile(toAppend):
            if "wgs" in currFile or "nucl" in currFile:
                nucl_accession_taxid_mapping_files.append(toAppend)
            elif "prot" in currFile:
                prot_accession_taxid_mapping_files.append(toAppend)

    # create 1 big file, for the reads and contigs mapped to their taxid, based on their accession_num
    start = time.time()
    contigs_reads_taxids_temp = analysis_path + "/contigs_reads_accessions_taxids_temp.txt"
    if os.path.isfile(contigs_reads_taxids_temp):
        os.remove(contigs_reads_taxids_temp)
    if join_seq_to_taxid(best_nt_scores, best_nr_scores, nucl_accession_taxid_mapping_files, prot_accession_taxid_mapping_files, contigs_reads_taxids_temp, analysis_path) is False: return None

    contigs_reads_taxids = analysis_path + "/contigs_reads_accessions_taxids.txt"
    subprocess.run("sed 's/ /\t/g' " + contigs_reads_taxids_temp + " | LC_COLLATE=C sort -k 1  > " + contigs_reads_taxids, shell=True) # change space to tabs
    os.remove(contigs_reads_taxids_temp)
    end = time.time()
    print("taxid identification via accessions took: " + str(end - start))

    # map reads to their contigs
    reads_belonging_to_contigs_unsorted = analysis_path + "/reads_belonging_to_contigs_unsorted.txt"
    get_best_sam_hits(reads_mapped_to_contigs_file, reads_belonging_to_contigs_unsorted, analysis_path)

    # now lets map the contig taxids to their reads
    # firstly retrieve all the dna/protein aligned reads

    species_avg_alignment_scores = analysis_path + "/species_avg_alignment_scores.txt"
    combined_nt_nr_scores_unsorted = analysis_path + "/combined_nt_nr_scores_unsorted.txt"
    subprocess.run("cat " + best_nt_scores + " > " + combined_nt_nr_scores_unsorted, shell=True)
    subprocess.run("cat " + best_nr_scores + " >> " + combined_nt_nr_scores_unsorted, shell=True)

    combined_nt_nr_scores = analysis_path + "/combined_nt_nr_scores.txt"
    subprocess.run("LC_COLLATE=C sort -k 1 " + combined_nt_nr_scores_unsorted + " > " + combined_nt_nr_scores, shell=True)
    match_scores(contigs_reads_taxids, combined_nt_nr_scores, analysis_path, species_avg_alignment_scores, args.taxdump_location)
    os.remove(combined_nt_nr_scores_unsorted)

    reads_belonging_to_contigs = analysis_path + "/reads_belonging_to_contigs.txt"
    subprocess.run("LC_COLLATE=C sort -k 2 " + reads_belonging_to_contigs_unsorted + " > " + reads_belonging_to_contigs, shell=True)
    os.remove(reads_belonging_to_contigs_unsorted)

    # stores reads that assembled, and its contigs taxid
    assembled_reads_taxids_temp = analysis_path + "/assembled_reads_taxids_temp.txt"
    contig_taxid_assignments = fetch_contig_taxids(contigs_reads_taxids)
    assign_taxids_to_assembled_reads(contig_taxid_assignments, reads_belonging_to_contigs, assembled_reads_taxids_temp)
    all_reads_taxids = analysis_path + "/all_reads_taxids.txt"
    subprocess.run("cat " + assembled_reads_taxids_temp + " > " + all_reads_taxids, shell=True) # change space to tabs

    # we need to now combine the contig aligned reads to the non contig aligned reads
    subprocess.run("awk '$1 !~ /^k[0-9]*/ {print $1\"\t\"$3}' " + contigs_reads_taxids + " >> " + all_reads_taxids, shell=True)

    # before proceeeding any further, lets remove the contaminants, simply remove the taxids
    contaminant_removal = True
    decontam_dirpath = ""

    if args.decontamination_input_files is None:
        print("Path to kraken index is not provided, skipping contamination removal")
        contaminant_removal = False
    else:
        if os.path.isdir(args.decontamination_input_files) is False:
            print("Directory that stores contaminant sequences does not exist, skipping contamination removal")
            contaminant_removal = False
        else:  
            decontam_dirpath = args.decontamination_input_files

            if decontam_dirpath[-1] == "/":
                decontam_dirpath = decontam_dirpath[0:-1]

            if (os.path.isdir(decontam_dirpath + "/others") is False and os.path.isdir(decontam_dirpath + "/controls") is False) or (os.path.listdir(decontam_dirpath + "/others") == 0 and os.path.listdir(decontam_dirpath + "/controls") == 0):
                print("No kraken outputs detected in " + decontam_dirpath + ". Please read manual regarding --decontamination_input_files. Skipping contamination removal step")
                contaminant_removal = False

    contaminants = []
    rcf_out = analysis_path + "/rcf_out.txt"
    if contaminant_removal is True:
        print("Starting decontamination")
        others = args.other_sequences
        others.append(args.inp1)
        others.append(args.inp2)
        start = time.time()
        contaminants = remove_contaminants_control(args.control_sequences, args.other_sequences, rcf_out, args.taxdump_location, analysis_path)
        end = time.time()
        print("Contaminant identification took: " + str(end - start))

    # now we can calculate read count method
    readCountsOutfile = analysis_path + "/readCountsOut.txt"

    # now lets get plot the abundances as krona charts
    # need to find lineages first
    all_taxids = fetch_taxids(contigs_reads_taxids)
    taxid_lineages_resolved, bad_taxids = get_lineage_info(all_taxids, args.taxdump_location)

    # bad_taxids are taxids whose ranks are not species (its higher than a species)

    countReads(all_reads_taxids, num_reads, readCountsOutfile, contaminants, bad_taxids)

    abundancesReadMethod = analysis_path + "/abundancesReadMethod.txt"
    wf = open(abundancesReadMethod, "w")
    with open(readCountsOutfile, "r") as f:
        for line in f:
            curr = line.split("\t")
            taxid = curr[0]
            score = curr[1].strip()
            if str(taxid) in taxid_lineages_resolved.keys():
                wf.write(score + "\t" + "\t".join(taxid_lineages_resolved[str(taxid)]) + "\n")

    wf.close()

    abundancesKronaReadMethod = final_plots_path + "/abundancesKronaReadMethod.html"
    subprocess.run("ImportText.pl " + abundancesReadMethod + " -o " + abundancesKronaReadMethod + " -fil " + species_avg_alignment_scores, shell=True)

    # now lets do abundance calculations via the tpm method
    # firstly get the length of the contigs
    # then get the number of reads aligned to them

    # firstly lets get the number of reads assigned to each contig
    contig_counts = getContigReadCount(reads_belonging_to_contigs)

    contig_unaligned_read_counts_temp = analysis_path + "/contig_unaligned_read_counts_temp.txt"
    if os.path.isfile(contig_unaligned_read_counts_temp):
        os.remove(contig_unaligned_read_counts_temp)
    getReadsLength(contigs_fq, contig_unaligned_read_counts_temp, contig_counts, False)


    # then repeat process of unaligned reads
    unassembled_reads_fwd = dirpath + "/preprocessing/unassembled_reads_fwd.fq"
    getReadsLength(unassembled_reads_fwd, contig_unaligned_read_counts_temp, None, True)

    # lets join contig_unaligned_read_counts with their taxid
    contig_reads_taxid_assignments = fetch_contig_read_taxids(contigs_reads_taxids)
     

    contig_unaligned_read_counts_len_taxid = analysis_path + "/contig_unaligned_read_counts_len_taxid.txt"
    compile_tpm_input_info(contig_reads_taxid_assignments, contig_unaligned_read_counts_temp, contig_unaligned_read_counts_len_taxid)

    log_path = dirpath + "/megahit_out/log"
    # retreiving n50 score
    n50 = float(subprocess.check_output('tail -n 2 ' + log_path + ' | head -n 1', shell=True).decode('utf-8').strip('\n').split(' ')[-2:-1][0])

    # now we can finally calculate TPM/FPKM
    tpm_abundance_file = analysis_path + "/tpm_fpkm.txt"
    get_abundance(contig_unaligned_read_counts_len_taxid, num_reads, n50, tpm_abundance_file, contaminants, bad_taxids)

    tpmAbundances = analysis_path + "/tpmAbundances.txt"
    wf = open(tpmAbundances, "w")
    with open(tpm_abundance_file, "r") as f:
        for line in f:
            curr = line.split("\t")
            taxid = curr[0]
            score = curr[1].strip()
            if str(taxid) in taxid_lineages_resolved.keys():
                wf.write(score + "\t" + "\t".join(taxid_lineages_resolved[str(taxid)]) + "\n")
    
    wf.close()

    tpmAbundancesKrona = final_plots_path + "/tpmAbundancesKrona.html"
    subprocess.run("ImportText.pl " + tpmAbundances + " -o " + tpmAbundancesKrona  + " -fil " + species_avg_alignment_scores, shell=True)

    # need to now go through each important output file to get relevant statistics for the output html

    summaryFile = dirpath + "/summary/summary.txt"
    summaryFileWriter = open(summaryFile, "a")

    # get number of assembled reads, put it into the megahit box
    numAssembledReads = countNumLines(reads_belonging_to_contigs)
    summaryFileWriter("numAssembledReads\t" + str(numAssembledReads) + "\n")

    # put it in the alignments box
    numUniqueNTHits = subprocess.check_output('sort -k 1 ' + nt_alignments_file + ' | cut -f 1 | uniq | wc -l | cut -d \' \' -f 1', shell=True)
    summaryFileWriter("numUniqueNTHits\t" + str(numUniqueNTHits) + "\n")

    numUniqueNRHits = subprocess.check_output('sort -k 1 ' + nr_alignments_file + ' | cut -f 1 | uniq | wc -l | cut -d \' \' -f 1', shell=True)
    summaryFileWriter("numUniqueNRHits\t" + str(numUniqueNRHits) + "\n")

    # put it in the identify best hits box
    numBestNTHits = countNumLines(best_nt_scores)
    summaryFileWriter("numBestNTHits\t" + str(numBestNTHits) + "\n")

    numBestNRHits = countNumLines(best_nr_scores)
    summaryFileWriter("numBestNRHits\t" + str(numBestNRHits) + "\n")

    totalUniqueAccessions = int(subprocess.check_output('cat ' + best_nt_scores + ' ' + best_nr_scores + ' | sort -k 2 | cut -f 2 | uniq | wc -l | cut -d \' \' -f 1', shell=True))
    summaryFileWriter("totalUniqueAccessions\t" + str(totalUniqueAccessions) + "\n")

    # put it in the obtain taxonomies box
    numMappedAccessions = subprocess.check_output('sort -k 2 '+ contigs_reads_taxids + ' | cut -f 2 | uniq | wc -l | cut -d \' \' -f 1', shell=True)
    summaryFileWriter("numMappedAccessions\t" + str(numMappedAccessions) + "\n")

    numUniqueTaxids = subprocess.check_output('sort -k 3 '+ contigs_reads_taxids + ' | cut -f 3 | uniq | wc -l | cut -d \' \' -f 1', shell=True)
    summaryFileWriter("numUniqueTaxids\t" + str(numUniqueTaxids) + "\n")

    numContigsWithTaxids = subprocess.check_output("egrep '^k[0-9] " + contigs_reads_taxids + ' | wc -l | cut -d \' \' -f 1', shell=True)
    summaryFileWriter("numContigsWithTaxids\t" + str(numContigsWithTaxids) + "\n")

    numUnassembledReadsWithTaxids = subprocess.check_output("egrep -v '^k[0-9] " + contigs_reads_taxids + ' | wc -l | cut -d \' \' -f 1', shell=True)
    summaryFileWriter("numUnassembledReadsWithTaxids\t" + str(numUnassembledReadsWithTaxids) + "\n")

    numAssembledReadsWithTaxids = countNumLines(assembled_reads_taxids_temp)
    summaryFileWriter("numAssembledReadsWithTaxids\t" + str(numAssembledReadsWithTaxids) + "\n")

    # nothing worthwhile to put in get median scores

    # put in get lineages for each taxid
    summaryFileWriter("numTaxidsWithRankSpecies\t" + str(len(taxid_lineages_resolved)) + "\n")
    summaryFileWriter("numTaxidsWithoutRankSpecies\t" + str(len(bad_taxids)) + "\n")

    # get numAssembledReadsWith a taxid
    # get number of contigs without a taxid
    # maybe have an option to expand this?

    summaryFileWriter.close()