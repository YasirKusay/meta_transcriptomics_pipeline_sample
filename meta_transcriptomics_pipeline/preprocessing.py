import argparse
import subprocess
import operator
import time
import os
from meta_transcriptomics_pipeline.helpers import check_fail
from meta_transcriptomics_pipeline.separate_reads_by_size import separate_reads_by_size
from meta_transcriptomics_pipeline.get_lineage_info import get_lineage_info
from meta_transcriptomics_pipeline.count_num_seqs import countNumSeqs

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

def preprocessing(args: argparse.Namespace):
    dirpath = args.dirpath

    if dirpath[-1] == "/":
        dirpath = dirpath[0:-1]

    if os.path.exists(dirpath + "/preprocessing") is False:
        os.mkdir(dirpath + "/preprocessing")

    # making analysis path to store the readCounts after each step, along with
    # any other relevant info

    stepsSummary = dirpath + "/analysis/stepsSummary.txt"

    if os.path.exists(dirpath + "/analysis") is False:
        os.mkdir(dirpath + "/analysis")

    dirpath = dirpath + "/preprocessing"

    ##################### FASTP ########################

    qc1 = dirpath + "/fastp_1.fastq"
    qc2 = dirpath + "/fastp_2.fastq"

    fastp_command = "fastp" +\
                    " --in1 " + args.inp1 +\
                    " --in2 " + args.inp2 +\
                    " --out1 " + qc1 +\
                    " --out2 " + qc2 +\
                    " -b 100 -B 100 " +\
                    " --qualified_quality_phred  " + args.qualified_quality_phred +\
                    " --unqualified_percent_limit " + args.unqualified_percent_limit +\
                    " --length_required " + args.length_required +\
                    " --low_complexity_filter " +\
                    " --detect_adapter_for_pe" +\
                    " --thread " + str(args.threads)
                    # need to consider adapters, should we give the user a chance to add adatpers?
                    # -b -B, means we want our reads/pairs to be at most 100 bases

    start = time.time()
    new_command = subprocess.run(fastp_command, shell=True)
    if check_fail("fastp", new_command) is True: return None
    end = time.time()
    print("fastp took: " + str(end - start))

    #################################### STAR HOST ###########################

    star_prefix = dirpath + "/star_host_"

    star_command = "STAR --genomeDir " + args.star_host_index + " --runThreadN " + str(args.threads) +\
                    " --readFilesIn " + qc1 + " " + qc2 + " --outFileNamePrefix " + star_prefix +\
                    " --outFilterMultimapNmax 99999 --outFilterScoreMinOverLread 0.5 --outFilterMatchNminOverLread 0.5" +\
                    " --outFilterMismatchNmax 999 --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx" +\
                    " --outSAMattributes Standard --quantMode TranscriptomeSAM GeneCounts --clip3pNbases 0"
    
    new_command = subprocess.run(star_command, shell=True)
    if check_fail("STAR", new_command) is True: return None

    star_host_dedup1 = star_prefix + "Unmapped_1.fastq"
    star_host_dedup2 = star_prefix + "Unmapped_2.fastq"

    os.rename(star_prefix + "Unmapped.out.mate1", star_host_dedup1)
    os.rename(star_prefix + "Unmapped.out.mate2", star_host_dedup2)

    #################################### SNAP HOST ###########################
    snap_host_mapping = dirpath + "/snap_host_mapped.bam"
    snap_host_command = "snap-aligner" + " paired " + args.snap_host_index + " " + star_host_dedup1 + " " + star_host_dedup2 +\
                    " -o " + snap_host_mapping + " -t " + str(args.threads) + " -I "
    start = time.time()
    new_command = subprocess.run(snap_host_command, shell=True)
    if check_fail("snap-aligner", new_command) is True: return None
    end = time.time()
    print("host subtraction via snap took: " + str(end - start))

    host_subtract_1 = dirpath + "/host_depleted1.fastq"
    host_subtract_2 = dirpath + "/host_depleted2.fastq"
    host_spare = dirpath + "/undetermined.fastq"

    # retrieving only host reads
    # using flag 12.
    # this flag means that we want reads where both its first and second pair failed to map
    # therefore, this treats reads pairs where only 1 read maps as host reads
    samtools_host_subtract_command = "samtools" + " fastq  -f 12 -@ " + str(args.threads) +\
                        " -1 " + host_subtract_1 +\
                        " -2 " + host_subtract_2 +\
                        " -s " + host_spare + " " +\
                        snap_host_mapping

    new_command = subprocess.run(samtools_host_subtract_command, shell=True)
    if check_fail("samtools", new_command) is True: return None

    #################################### SORTMERNA ############################
    aligned = dirpath + "/aligned"
    noRna = dirpath + "/noRrna"
    noRna1 = dirpath + "/noRrna_fwd.fq"
    noRna2 = dirpath + "/noRrna_rev.fq"

    sortmerna_command = "sortmerna" +\
                    " --ref " + args.sortmerna_rrna_database +\
                    " --aligned " + aligned +\
                    " --other " + noRna +\
                    " --fastx " +\
                    " --reads " + host_subtract_1 + " --reads " + host_subtract_2 +\
                    " --threads " + str(args.threads) +\
		            " --out2 TRUE " +\
		            " --paired_in TRUE"

    start = time.time()
    new_command = subprocess.run(sortmerna_command, shell=True)
    if check_fail("sortmerna", new_command) is True: return None
    end = time.time()
    print("sortmerna took: " + str(end - start))

    #################################### CLUMPIFY DEDUP #######################

    fullyQc1 = dirpath + "/fullyQc_fwd.fq"
    fullyQc2 = dirpath + "/fullyQc_rev.fq"

    clumpify_command = "clumpify.sh  in1=" + noRna1 + " in2=" + noRna2 +\
                            " out1=" + fullyQc1 + " out2=" + fullyQc2 + " dedupe=t"
    
    new_command = subprocess.run(clumpify_command, shell=True)
    if check_fail("clumpify.sh", new_command) is True: return None

    # QUICK ALIGNMENT, JUST ALIGN REMAINING READS USING KRAKEN AGAINST KRAKEN_PLUS
    '''
    num_reads_bytes = subprocess.run(['grep', '-c', '.*', fullyQc1], stdout=subprocess.PIPE)
    num_reads_str = num_reads_bytes.stdout.decode('utf-8')
    num_reads = int(num_reads_str.replace('\n', ''))/4 # finally in int format, dividing by 4 because its in fastq format
    fast_mode_output = dirpath + "/fast_mode_output"
    kraken_command = "kraken2 --db " + args.kraken_db + " --threads " + str(args.threads) +\
                        " --output " + fast_mode_output + " --paired " + fullyQc1 + " " + fullyQc2
    
    new_command = subprocess.run(kraken_command, shell=True)    
    if check_fail("kraken", new_command, []) is True: return None
    kraken_res_out = dirpath + "/kraken_res_out"
    process_fast_mode_output(fast_mode_output, kraken_res_out, num_reads)
    fastAbundances = dirpath + "/fastAbundances.txt"
    get_lineage_info(kraken_res_out, fastAbundances, args.taxdump_location)
    fastAbundancesKrona = dirpath + "/fastAbundancesKrona.html"
    subprocess.run("ImportText.pl " + fastAbundances + " -o " + fastAbundancesKrona, shell=True)
    '''
    
    #################### MEGAHIT ###########################
    megahit_out_path = args.dirpath

    if megahit_out_path[-1] == "/":
        megahit_out_path = megahit_out_path[0:-1]
    megahit_out_path = megahit_out_path + "/megahit_out"

    megahit_command = "megahit" + " -1 " + fullyQc1 +\
                        " -2 " + fullyQc2 +\
                        " -o " + megahit_out_path + " -t " + str(args.threads) # is an output directory 
    contigs = megahit_out_path + "/final_contigs.fa"
    start = time.time()
    new_command = subprocess.run(megahit_command, shell=True)
    if check_fail("megahit", new_command) is True: return None
    end = time.time()
    print("assembly via megahit took: " + str(end - start))
    new_command = subprocess.run("mv " + megahit_out_path + "/final.contigs.fa " + contigs, shell=True)

    # we must retrieve the unaligned reads
    reads_mapped_to_contigs_file = dirpath + "/reads_mapped_to_contigs.sam"
    align_reads_to_contigs_cmd = "bbwrap.sh" + " ref=" + contigs +\
                                " in=" + fullyQc1 +\
                                " in2=" + fullyQc2 +\
                                " -out=" + reads_mapped_to_contigs_file  
    new_command = subprocess.run(align_reads_to_contigs_cmd, shell=True)
    if check_fail("megahit", new_command) is True: return None 

    # now lets retrieve the reads that did not align
    unassembled_reads_fwd = dirpath + "/unassembled_reads_fwd.fq"
    unassembled_reads_rev = dirpath + "/unassembled_reads_rev.fq"

    # same principle here as the host mapping step
    align_command = "samtools fastq -f 12 -1 " + unassembled_reads_fwd +\
                    " -2 " + unassembled_reads_rev + " " + reads_mapped_to_contigs_file
    new_command = subprocess.run(align_command, shell=True)
    if check_fail("samtools", new_command) is True: return None 

    unassembled_reads_shorter_1 = dirpath + "/unassembled_reads_shorter_fwd.fq"
    unassembled_reads_shorter_2 = dirpath + "/unassembled_reads_shorter_rev.fq"
    unassembled_reads_longer_1 = dirpath + "/unassembled_reads_longer_fwd.fq"
    unassembled_reads_longer_2 = dirpath + "/unassembled_reads_longer_rev.fq"

    separate_reads_by_size(unassembled_reads_fwd, unassembled_reads_longer_1, unassembled_reads_shorter_1)
    separate_reads_by_size(unassembled_reads_rev, unassembled_reads_longer_2, unassembled_reads_shorter_2)

    # need to convert file above from fa to fq, simply done using seqtk
    contigs_fq = megahit_out_path + "/final_contigs.fq"
    seqtk_command = "seqtk" + " seq -F '#' " + contigs + " > " + contigs_fq
    new_command = subprocess.run(seqtk_command, shell=True)
    if check_fail("seqtk", new_command) is True: return None

    # need to merge paired end reads
    merged_pe = dirpath + "/merged_reads.fq"
    merge_command = "seqtk mergepe " + unassembled_reads_fwd + " " + unassembled_reads_rev + " > " + merged_pe
    new_command = subprocess.run(merge_command, shell=True)
    if check_fail("seqtk mergepe", new_command) is True: return False 

    # merge pe reads with contigs
    combined_file_fq = dirpath + "/combined_reads_contigs_file.fq"
    subprocess.run("cat " + contigs_fq + " > " + combined_file_fq, shell=True)
    subprocess.run("cat " + merged_pe + " >> " + combined_file_fq, shell=True)

    # need to convert above to fasta
    combined_file_fa = dirpath + "/combined_reads_contigs_file.fa"
    seqtk_command = "seqtk" + " seq -a " + combined_file_fq + " > " + combined_file_fa
    new_command = subprocess.run(seqtk_command, shell=True)
    if check_fail("seqtk", new_command) is True: return False 

    # creating alignment folder
    alignments_path = args.dirpath

    if alignments_path[-1] == "/":
        alignments_path = alignments_path[0:-1]
    
    alignments_path = alignments_path + "/alignments"

    if os.path.exists(alignments_path) is False:
        os.mkdir(alignments_path)

    # need to now go through each important output file to get relevant statistics for the output html

    summaryPath = dirpath + "/summary"
    if os.path.exists(summaryPath) is False:
        os.mkdir(summaryPath)

    summaryFile = summaryPath + "/summary.txt"
    summaryFileWriter = open(summaryFile, "w")

    numReadsAtStart = 0
    if (os.path.splitext(args.inp1) == ".fa" or os.path.splitext(args.inp1) == ".fna" or os.path.splitext(args.inp1) == ".fasta"):
        numReadsAtStart =  countNumSeqs(args.inp1, True)
    else:
        numReadsAtStart =  countNumSeqs(args.inp1, False)

    summaryFileWriter("Start\t" + str(numReadsAtStart) + "\n")
    
    numReadsAfterFastq = countNumSeqs(qc1, False)
    summaryFileWriter("Fastq\t" + str(numReadsAfterFastq) + "\n")
    
    numReadsAfterStar = countNumSeqs(star_host_dedup1, False)
    summaryFileWriter("Star\t" + str(numReadsAfterStar) + "\n")
    
    numReadsAfterSnap = countNumSeqs(host_subtract_1, False)
    summaryFileWriter("Snap\t" + str(numReadsAfterSnap) + "\n")
    
    numReadsAfterSortmerna = countNumSeqs(noRna1, False)
    summaryFileWriter("Sortmerna\t" + str(numReadsAfterSortmerna) + "\n")
    
    numReadsAfterClumpify = countNumSeqs(fullyQc1, False)
    summaryFileWriter("Clumpify\t" + str(numReadsAfterClumpify) + "\n")

    log_path = megahit_out_path + "/log"
    assemblyStats = subprocess.check_output('tail -n 2 ' + log_path + ' | head -n 1', shell=True).decode('utf-8').strip('\n').split(' ')

    numContigs = assemblyStats[2]
    summaryFileWriter("numContigs\t" + str(numContigs) + "\n")

    totalBases = assemblyStats[5]
    summaryFileWriter("totalContigsBases\t" + str(totalBases) + "\n")

    shortestContig = assemblyStats[8]
    summaryFileWriter("shortestContig\t" + str(shortestContig) + "\n")
    
    longestContig = assemblyStats[11]
    summaryFileWriter("longestContig\t" + str(longestContig) + "\n")

    avgContigLength = assemblyStats[14]
    summaryFileWriter("avgContigLength\t" + str(avgContigLength) + "\n")

    n50 = assemblyStats[17]
    summaryFileWriter("n50\t" + str(n50) + "\n")

    # will retreive numAssembledContigs at the finalisation step
    numReadsUnassembled = countNumSeqs(unassembled_reads_fwd, False)
    summaryFileWriter("numReadsUnassembled\t" + str(numReadsUnassembled) + "\n")

    numLongReadsUnassembled = countNumSeqs(unassembled_reads_longer_1, False)
    summaryFileWriter("numLongReadsUnassembled\t" + str(numLongReadsUnassembled) + "\n")

    numShortReadsUnassembled = countNumSeqs(unassembled_reads_shorter_1, False)
    summaryFileWriter("numShortReadsUnassembled\t" + str(numShortReadsUnassembled) + "\n")

    summaryFileWriter.close()