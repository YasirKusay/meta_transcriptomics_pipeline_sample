import argparse
import subprocess
import os
import time
#from main_pipeline.helpers import check_command_exists, check_fail, generate_temp_file
from helpers import check_command_exists, check_fail, generate_temp_file
from merge_sams import merge_sams
from merge_contigs import merge_contigs
from join_taxid_contigs import join_taxid_contigs

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

def re_adjust_tpms(file, numlines, path):
    total_lines = numlines
    total_original = 0

    new_tpm_file = path + "/new_tpm_file"

    with open(file, "r") as f:
        for line in f:
            curr = line.split()
            if (float(curr[1]) == 0):
                total_lines = total_lines + 1
            else:
                total_original = total_original + float(curr[1])

    to_subtract = total_lines/total_original
    wf = open(new_tpm_file, "w")

    with open(file, "r") as f:
        for line in f:
            curr = line.split()
            if (float(curr[1]) == 0):
                curr[1] = "1"
            else:
                curr[1] = str(float(curr[1]) - to_subtract)

            wf.write(curr[0] + "\t" + curr[1] + "\n")

    wf.close()
    return new_tpm_file

def get_abundance(joined, final_file):
    wf = open(final_file, "w")
    final_res = {}
    with open(joined, "r") as f:
        for line in f:
            curr = line.split()
            if curr[3] in final_res.keys():
                final_res[curr[3]] += float(curr[2])
            else:
                final_res[curr[3]] = float(curr[2])

    for key in final_res.keys():
        wf.write(key + "\t" + str(final_res[key]) + "\n")

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

    out1 = generate_temp_file(".fastq", dirpath)
    out2 = generate_temp_file(".fastq", dirpath)
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
                    " --thread " + str(args.threads)
                    # need to consider adapters, should we give the user a chance to add adatpers?

    start = time.time()
    new_command = subprocess.run(fastp_command, shell=True)
    if check_fail(fastp_path, new_command, [out1, out2]) is True: return None
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
    new_command = subprocess.run(sortmerna_command, shell=True)
    if check_fail(sortmerna_path, new_command, [out1, out2, aligned + "_fwd.fastq", aligned + "_rev.fastq",  other + "_fwd.fastq", other + "_rev.fastq"]) is True: return None
    end = time.time()
    print("sortmerna took: " + str(end - start))
    #os.remove(aligned + "_fwd.fastq", aligned + ".log", aligned + "_rev.fastq")
    generated_files.append(other + "_fwd.fastq")
    generated_files.append(other + "_rev.fastq")

    ############################# SNAP HUMAN #############################
    human_out = dirpath + "/snap_human_out.sam"
    snap_path = 'snap-aligner'
    snap_human_command = snap_path + " paired " + args.snap_human_index + " " + other + "_fwd.fastq " + other + "_rev.fastq " +\
                    " -o " + human_out + " -t " + str(args.threads)
    start = time.time()
    new_command = subprocess.run(snap_human_command, shell=True)
    if check_fail(snap_path, new_command, [generated_files]) is True: return None
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

    new_command = subprocess.run(samtools_human_subtract_command, shell=True)
    if check_fail(samtools_path, new_command, []) is True: return None
    generated_files.append(human_subtract_1 + ".fastq")
    generated_files.append(human_subtract_2 + ".fastq")
    generated_files.append(human_spare + ".fastq")

    #################### MEGAHIT ###########################
    megahit_path = "megahit"
    contig_path = dirpath + "/megahit_out"
    megahit_command = megahit_path + " -1 " + human_subtract_1 +\
                        " -2 " + human_subtract_2 +\
                        " -o " + contig_path + " -t " + str(args.threads) # is an output directory 
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
    #if run_snap_single(args.snap_index, new_contigs, snap_contigs, args.threads) == False: return None

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

    # we must retrieve the unaligned reads
    bbwrap_path = "bbwrap.sh"
    reads_mapped_to_contigs_file = dirpath + "/reads_mapped_to_contigs.sam"
    align_reads_to_contigs_cmd = bbwrap_path + " ref=" + contigs +\
                                " in=" + human_subtract_1 +\
                                " in2=" + human_subtract_2 +\
                                " -out=" + reads_mapped_to_contigs_file  
    new_command = subprocess.run(align_reads_to_contigs_cmd, shell=True)
    if check_fail(bbwrap_path, new_command, []) is True: return None 

    # now lets retrieve the reads that did not align
    new_fwd = dirpath + "/new_fwd.fq"
    new_rev = dirpath + "/new_rev.fq"

    align_command = "samtools fastq -f4 -1 " + new_fwd +\
                    " -2 " + new_rev + " " + reads_mapped_to_contigs_file
    new_command = subprocess.run(align_command, shell=True)
    if check_fail(samtools_path, new_command, []) is True: return None 

    # running exact set of commands above but for unaligned reads this time

    snap_reads = dirpath + "/snap_reads_out.sam"
    snap_read_command = "snap-aligner" + " paired " + args.snap_index + " " + new_fwd + " " + new_rev +\
                    " -o " + snap_reads + " -t " + str(args.threads)
    start = time.time()
    new_command = subprocess.run(snap_read_command, shell=True)
    if check_fail("snap-aligner", new_command, []) is True: return False
    end = time.time()
    print("read alignment against nt took: " + str(end - start))

    # now lets align reads
    # need to merge paired end reads first though
    merged_pe = dirpath + "/merged_reads.fq"
    merge_command = "seqtk mergepe " + new_fwd + " " + new_rev + " > " + merged_pe
    new_command = subprocess.run(merge_command, shell=True)
    if check_fail("seqtk mergepe", new_command, []) is True: return False 

    diamond_reads = dirpath + "/diamond_reads_out.sam"
    start = time.time()
    if run_diamond(args.diamond_index, merged_pe, diamond_reads, args.threads, 101) == False: return None
    end = time.time()
    print("read alignment against nr took: " + str(end - start))

    

    snap_reads_out, diamond_reads_out = merge_sams(snap_reads, diamond_reads, dirpath)

    contigs_reads_taxids = dirpath + "/nucl_prot_taxids.txt"
    if join_taxid_contigs(snap_contig_out, snap_reads_out, diamond_contig_out, diamond_reads_out, args.nucl_accession_taxid_mapping, args.prot_accession_taxid_mapping, contigs_reads_taxids, dirpath) is False: return None

    # running tpm calculations via rsem
    # need to firstly build an index
    path_to_rsem_index = dirpath + "/rsem_ref"
    rsem_command_prep = "rsem-prepare-reference --bowtie2 " + " " + contigs + " " + path_to_rsem_index
    new_command = subprocess.run(rsem_command_prep, shell=True)
    if check_fail("rsem-prepare-reference", new_command, []) is True: return False 

    # now lets run rsem
    rsem_out = dirpath + "/rsem_out"
    rsem_command_calc = "rsem-calculate-expression --paired-end --bowtie2 " +\
                        human_subtract_1 + " " + human_subtract_2 + " " +\
                        path_to_rsem_index + " " + rsem_out

    new_command = subprocess.run(rsem_command_calc, shell=True)
    if check_fail("rsem-calculate-expression", new_command, []) is True: return False 

    gene_file = rsem_out + ".genes.results"

    # now lets rescale the tpm's
    all_unmapped = dirpath + "/all_unmapped"
    cmd1 = 'egrep "@" ' + new_fwd + " > " + all_unmapped
    new_command = subprocess.run(cmd1, shell=True)
    cmd2 = 'egrep "@" ' + new_rev + " >> " + all_unmapped
    new_command = subprocess.run(cmd2, shell=True)
    all_sorted = dirpath + "/all_sorted"
    cmd3 = "sort -k1 " + all_unmapped + " | uniq | sed 's/@//' > " + all_sorted
    new_command = subprocess.run(cmd3, shell=True)
    all_sorted_tpm = dirpath + "/all_sorted_tpm" # tpm refers to transcripts per million
    add_p = "awk 'BEGIN { FS = OFS = \"\t\" } {print($0, \"1\")}' " + all_sorted + " > " + all_sorted_tpm
    new_command = subprocess.run(add_p, shell=True)
    num_lines = sum(1 for line in open(all_sorted))

    combined_tpms = dirpath + "/combined_tpms"
    new_command = subprocess.run("awk 'NR > 1' " + gene_file + " | cut -f1,6 | sort -k2 -n > " + combined_tpms, shell=True)
    final_tpms = re_adjust_tpms(combined_tpms, num_lines, dirpath)

    new_command = subprocess.run("cat " + all_sorted_tpm + " >> " + final_tpms, shell=True)
    final_tpms_sorted = dirpath + "/final_tpms_sorted"
    new_command = subprocess.run("sort -k1 " + final_tpms + " > " + final_tpms_sorted, shell=True)
    contigs_reads_taxids_sorted = dirpath + "/nucl_prot_taxids_sorted.txt"
    new_command = subprocess.run("sort -k1 " + contigs_reads_taxids + " > " + contigs_reads_taxids_sorted, shell=True)

    join_file = dirpath + "/join_file"
    join_tpms_taxids = "join -1 1 -2 1 -o \"2.1 2.2 1.2 2.3\" " + final_tpms_sorted + " " + contigs_reads_taxids_sorted + " > " + join_file
    new_command = subprocess.run(join_tpms_taxids, shell=True)

    final_res = dirpath + "/final_res" 
    get_abundance(join_file, final_res)

    print("DONE!!!")

    # we are done, lets remove the temp directory
    #shutil.rmtree(dirpath)

    
    print("DONE!!!")