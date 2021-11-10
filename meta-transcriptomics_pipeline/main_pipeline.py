import argparse
import subprocess
import os
#from main_pipeline.helpers import check_command_exists, check_fail, generate_temp_file
from helpers import check_command_exists, check_fail, generate_temp_file
from merge_sams import merge_sams
from join_taxid_contigs import join_taxid_contigs

def run_snap_single(index, in_path, out_path, threads):
    snap_contig_command = "snap-aligner" + " single " + index + " " + in_path +\
                    " -o " + out_path + " -t " + str(threads)
    new_command = subprocess.run(snap_contig_command, shell=True)
    if check_fail("snap-aligner", new_command, []) is True: return False
    
    return True

def run_diamond(index, in_path, out_path, threads):
    diamond_command = "diamond" + " blastx -db " + index +\
                        " --query " + in_path + " --sensitive --max-target-seqs 1 --outfmt 101" +\
                        " --threads " + str(threads) +\
                        " --out " + out_path
    new_command = subprocess.run(diamond_command, shell=True)
    if check_fail("diamond", new_command, []) is True: return False

    return True


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

    new_command = subprocess.run(fastp_command, shell=True)
    if check_fail(fastp_path, new_command, [out1, out2]) is True: return None

    print("PART 1 DONE FASTP")
    print(out1)
    print(out2)
    print(os.path.isfile(out1))
    print(new_command.returncode)
    

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

    print("PART 2 SORTME")
    print(aligned)
    print(other)
    print(os.path.isfile(out1))
    print(os.path.isfile(out2))
    print(args.sortmerna_index)

    new_command = subprocess.run(sortmerna_command, shell=True)
    if check_fail(sortmerna_path, new_command, [out1, out2, aligned + "_fwd.fastq", aligned + "_rev.fastq",  other + "_fwd.fastq", other + "_rev.fastq"]) is True: return None
    #os.remove(aligned + "_fwd.fastq", aligned + ".log", aligned + "_rev.fastq")
    generated_files.append(other + "_fwd.fastq")
    generated_files.append(other + "_rev.fastq")

    #if check_fail(new_command, [file1, file2]) is False: return None
    print("PART 2 DONE SORTME")
    print(new_command.returncode)
    print(os.path.isfile(out1))

    print('out: ', new_command.stdout)
    print('err: ', new_command.stderr)

    ############################# SNAP HUMAN #############################
    human_out = dirpath + "/snap_human_out.sam"
    snap_path = 'snap-aligner'
    snap_human_command = snap_path + " paired " + args.snap_human_index + " " + other + "_fwd.fastq " + other + "_rev.fastq " +\
                    " -o " + human_out + " -t " + str(args.threads)
    new_command = subprocess.run(snap_human_command, shell=True)
    if check_fail(snap_path, new_command, [generated_files]) is True: return None

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
    new_command = subprocess.run(megahit_command, shell=True)
    if check_fail(megahit_path, new_command, []) is True: return None
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
    snap_contigs = dirpath + "/snap_contigs_out.sam"
    if run_snap_single(args.snap_index, new_contigs, snap_contigs, args.threads) == False: return None

    # now lets align reads
    diamond_contigs = dirpath + "/diamond_contigs_out.sam"
    if run_diamond(args.diamond_index, new_contigs, diamond_contigs, args.threads) == False: return None

    snap_contig_out, diamond_contig_out = merge_sams(snap_contigs, diamond_contigs, dirpath)

    # we must retrieve the unaligned reads
    bbwrap_path = "bbwrap.sh"
    reads_mapped_to_contigs_file = dirpath + "/reads_mapped_to_contigs.sam"
    align_reads_to_contigs_cmd = bbwrap_path + " ref=" + contigs +\
                                " in=" + human_subtract_1 + ".fastq" +\
                                " in2=" + human_subtract_2 + ".fastq" +\
                                " -out=" + reads_mapped_to_contigs_file  
    new_command = subprocess.run(align_reads_to_contigs_cmd)
    if check_fail(bbwrap_path, new_command, []) is True: return None 

    # now lets retrieve the reads that did not align
    new_fwd = dirpath + "/new_fwd.fq"
    new_rev = dirpath + "/new_rev.fq"

    align_command = "samtools fastq -f4 -1 " + new_fwd +\
                    " -2 " + new_rev + " " + reads_mapped_to_contigs_file
    new_command = subprocess.run(align_command)
    if check_fail(samtools_path, new_command, []) is True: return None 

    # running exact set of commands above but for unaligned reads this time

    snap_reads = dirpath + "/snap_reads_out.sam"
    snap_read_command = "snap-aligner" + " paired " + args.snap_index + " " + new_fwd + " " + new_rev +\
                    " -o " + snap_reads + " -t " + str(args.threads)
    new_command = subprocess.run(snap_read_command, shell=True)
    if check_fail("snap-aligner", new_command, []) is True: return False

    # now lets align reads
    # need to merge paired end reads first though
    merged_pe = dirpath + "merged_reads.fq"
    merge_command = "seqtk mergepe " + new_fwd + " " + new_rev + " > " + merged_pe
    new_command = subprocess.run(merge_command, shell=True)
    if check_fail("seqtk mergepe", new_command, []) is True: return False 

    diamond_reads = dirpath + "/diamond_reads_out.sam"
    if run_diamond(args.diamond_index, merged_pe, diamond_reads, args.threads) == False: return None

    snap_reads_out, diamond_reads_out = merge_sams(snap_reads, diamond_reads, dirpath)

    contigs_reads_taxids = dirpath + "/nucl_prot_taxids.txt"
    if join_taxid_contigs(snap_contig_out, snap_reads_out, diamond_contig_out, diamond_reads_out, args.nucl_accession_taxid_mapping, args.prot_accession_taxid_mapping, contigs_reads_taxids, dirpath) is False: return None

    # running tpm calculations via rsem
    # need to firstly build an index


    print("DONE!!!")

    # we are done, lets remove the temp directory
    #shutil.rmtree(dirpath)

