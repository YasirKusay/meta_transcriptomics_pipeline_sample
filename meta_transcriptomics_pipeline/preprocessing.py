import argparse
import math
import subprocess
import operator
import time
import os
import json
import shutil
import numpy as np
import matplotlib.pyplot as plt
from meta_transcriptomics_pipeline.helpers import run_shell_command
from meta_transcriptomics_pipeline.separate_reads_by_size import separate_reads_by_size
from meta_transcriptomics_pipeline.get_lineage_info import get_lineage_info
from meta_transcriptomics_pipeline.count_num_lines import countNumSeqs

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

    if os.path.exists(dirpath + "/analysis") is False:
        os.mkdir(dirpath + "/analysis")

    if os.path.exists(dirpath + "/final_plots") is False:
        os.mkdir(dirpath + "/final_plots")

    dirpath = dirpath + "/preprocessing"

    # failedStep = None
    # failFile = dirpath + "/failedStep.txt"

    # # records the part that the preprocessing step failed at
    # # hence, we do not need to rerun the entire pipeline should just one step fail
    # if os.path.exists(failFile):
    #     if args.continueFromFailure is True and os.stat(failFile).st_size > 0:
    #         with open(failFile, "r") as f:
    #             failedStep = f.readline().strip()
    #     else:
    #         os.remove(failFile)

    ##################### FASTP ########################

    qc1 = dirpath + "/fastp_1.fastq"
    qc2 = dirpath + "/fastp_2.fastq"
    fastp_json = dirpath + "/fastp.json"

    fastp_command = "fastp" +\
                    " --in1 " + args.inp1 +\
                    " --in2 " + args.inp2 +\
                    " --out1 " + qc1 +\
                    " --out2 " + qc2 +\
                    " --json " + fastp_json +\
                    " -b " + args.read_length + " -B " + args.read_length +\
                    " --qualified_quality_phred  " + args.qualified_quality_phred +\
                    " --unqualified_percent_limit " + args.unqualified_percent_limit +\
                    " --length_required " + args.length_required +\
                    " --low_complexity_filter " +\
                    " --detect_adapter_for_pe" +\
                    " --thread " + str(args.threads)
                    # need to consider adapters, should we give the user a chance to add adatpers?
                    # -b -B, means we want our reads/pairs to be at most 100 bases

    start = time.time()
    run_shell_command(fastp_command)
    end = time.time()
    print("Fastp took: " + str(end - start))

    #################################### STAR HOST ###########################

    star_prefix = dirpath + "/star_host_"

    star_command = "STAR --genomeDir " + args.star_host_index + " --runThreadN " + str(args.threads) +\
                    " --readFilesIn " + qc1 + " " + qc2 + " --outFileNamePrefix " + star_prefix +\
                    " --outFilterMultimapNmax 99999 --outFilterScoreMinOverLread 0.5 --outFilterMatchNminOverLread 0.5" +\
                    " --outFilterMismatchNmax 999 --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx" +\
                    " --outSAMattributes Standard --quantMode TranscriptomeSAM GeneCounts --clip3pNbases 0"
    
    start = time.time()
    run_shell_command(star_command)
    end = time.time()
    print("Human depletion via star took: " + str(end - start))

    star_host_dedup1 = star_prefix + "Unmapped_1.fastq"
    star_host_dedup2 = star_prefix + "Unmapped_2.fastq"

    os.rename(star_prefix + "Unmapped.out.mate1", star_host_dedup1)
    os.rename(star_prefix + "Unmapped.out.mate2", star_host_dedup2)

    # need to extract ERCC coverages (if they exist)
    if os.path.exists(dirpath + "/star_coverage.txt"):
        os.remove(dirpath + "/star_coverage.txt")
    if os.path.exists(dirpath + "/ercc_coverage.txt"):
        os.remove(dirpath + "/ercc_coverage.txt")

    erccReadCounts = None

    # check if the output file is empty
    if args.ercc_expected_concentration is not None and os.path.exists(args.ercc_expected_concentration) is True and os.stat(args.ercc_expected_concentration).st_size > 0:
        start = time.time()
        run_shell_command("pileup.sh in=" + dirpath + "/star_host_Aligned.sortedByCoord.out.bam" + " out=" + dirpath + "/star_coverage.txt" + " -Xmx" + str(args.memory) + "g" + " secondary=false")
        end = time.time()

        print("pileup.sh on the star output took: " + str(end - start))

        # running this manually as if we fail to grep anything, it will return NOT 1 and hence will exit.
        subprocess.run("egrep -e \"^ERCC\" " + dirpath + "/star_coverage.txt" + " > " + dirpath + "/ercc_coverage.txt", shell=True)

        # dict stores ercc as keys, and its values is a list with the expected as first value and readcount (actual value) as the second value
        # the first value within these sublists will be the "x" within the plot and the secon value will be the "y"
        ercc_counts = {}

        if os.path.exists(dirpath + "/ercc_coverage.txt") is True and os.stat(dirpath + "/ercc_coverage.txt").st_size > 0:

            # https://www.biostars.org/p/110052/
            # need to get the ERCC couunts, firstly grab the ERCC sequences from the bam file, get instances where both reads
            # are mapped (I strongly believe that the unmapped files will also include sequences where only one pair is mapped)
            # then count up the sequences

            regions = []

            with open(dirpath + "/ercc_coverage.txt", "r") as f:
                for line in f:
                    curr = line.split("\t")
                    regions.append(curr[0])

            run_shell_command("samtools index " + dirpath + "/star_host_Aligned.sortedByCoord.out.bam")
            run_shell_command("samtools view -F 260 " + dirpath + "/star_host_Aligned.sortedByCoord.out.sam " + " ".join(regions) + " -o " + dirpath + "/star_host_Aligned.ERCC_only.out.sam")
            erccReadCounts = int(subprocess.check_output("cut -f 1 " + dirpath + "/star_host_Aligned.ERCC_only.out.sam" + " | sort -k 1 | uniq | wc -l", shell=True).decode('utf-8').strip())

            # sort the ercc_expected_concentration file just in case, such that our ercc_coverage.txt file can be compared simultaneously
            with open(args.ercc_expected_concentration) as f:
                for line in f:
                    if line == "\n":
                        continue
                    if line == "":
                        break
                    curr = line.strip().split("\t")
                    if len(curr) < 2:
                        continue
                    if float(curr[1]) == 0:
                        continue
                    ercc_counts[curr[0]] = [math.log(float(curr[1]), 10)]

            # now do the same but for the actual read counts
            with open(dirpath + "/ercc_coverage.txt") as f:
                for line in f:
                    if line == "\n":
                        continue
                    if line == "":
                        break
                    curr = line.strip().split("\t")
                    if len(curr) < 2:
                        continue
                    if float(curr[1]) == 0:
                        continue
                    if curr[0] in ercc_counts.keys():
                        ercc_counts[curr[0]].append(math.log(float(curr[1]), 10))

            # in case we missed anything
            # actually remove them, treat as 0 for read count

            keys_to_remove = []

            for key in ercc_counts.keys():    
                if len(ercc_counts[key]) == 1:
                    keys_to_remove.append(key)

            for key in keys_to_remove:
                del ercc_counts[key]

            #for key in ercc_counts.keys():
            #    if len(ercc_counts[key]) == 1:
            #        ercc_counts[key].append(0)
            # check if we have at least one value

            # check if we have at least one value
            if len(ercc_counts) > 0:
                data = np.asarray(list(ercc_counts.values()))
                x_values = data[:,0]
                y_values = data[:,1]
                #plt.plot(data[:,0], data[:, 1])
                a, b = np.polyfit(x_values, y_values, 1)
                plt.scatter(x_values, y_values, color='red', edgecolor='black', linewidth=0.4)

                best_fit = [a*x_val+b for x_val in x_values]
                # below wont show dots spaced evenly: https://stackoverflow.com/questions/29906892/pyplot-cannot-draw-dotted-line
                # plt.plot(x_values, a*x_values+b, color='black', linestyle=(0, (1, 10)), linewidth=0.8)

                # calculating r^2: https://towardsdatascience.com/r-squared-recipe-5814995fa39a
                # starting at step 2
                mean_y = y_values.mean()

                # step 3
                differences_regression_line = []
                for i in range(0, len(x_values)):
                    differences_regression_line.append((a * x_values[i] + b) - y_values[i])

                regression_sum = 0
                for i in differences_regression_line:
                    regression_sum += (i*i)

                # step 4
                differences_mean_line = [mean_y - y_val for y_val in y_values]
                mean_sum = 0
                for i in differences_mean_line:
                    mean_sum += (i*i)

                # step 5
                r_squared = (mean_sum - regression_sum)/mean_sum

                plt.plot([min(x_values), max(x_values)], [min(best_fit), max(best_fit)], color='black', linestyle=':', linewidth=0.8)
                plt.text(min(x_values), max(y_values), "R² = " + str(r_squared))

                # plt.plot(x_values, a*x_values+b, label = "R² = ",color='black', linestyle='--', linewidth=2)
                # plt.text(, "R² = ") # trendline equation
                # we have got the above from: https://www.statology.org/line-of-best-fit-python/
                plt.xlabel("Log10 ERCC spike-in concentration")
                plt.ylabel("Log10 coverage per gene")
                plt.savefig(dirpath + '/../final_plots/ercc_plot.png', dpi=300) # dpi to control resolution
        else:
            print("The coverage file generated from the snap output (preprocessing/star_coverage.txt) is empty. ERCC step has failed but the pipeline will continue. In the future, please recheck the STAR index provided.")    
    else:
        print("Detected ERCC sequences in the star index but the file that will be used to compare the expected ercc concentration has not been provided/or is empty. Please provide this file to --ercc_expected_concentration.")

    #################################### SNAP HOST ###########################
    snap_host_mapping_unsorted = dirpath + "/snap_host_mapped_unsorted.bam"
    snap_host_command = "snap-aligner" + " paired " + args.snap_host_index + " " + star_host_dedup1 + " " + star_host_dedup2 +\
                    " -o " + snap_host_mapping_unsorted + " -t " + str(args.threads) + " -I "
    start = time.time()
    run_shell_command(snap_host_command)
    end = time.time()
    print("Host subtraction via snap took: " + str(end - start))

    snap_host_mapping = dirpath + "/snap_host_mapped.bam"

    # the snap output file must be sorted by name otherwise samtools fastq will not work as intended
    # biostars.org/p/303292/
    samtools_sort_command = "samtools sort -@ " + str(args.threads) + " -n " + snap_host_mapping_unsorted + " -o " + snap_host_mapping
    start = time.time()
    run_shell_command(samtools_sort_command)
    end = time.time()
    print("Sorting the sam file generated by the snap depletion step took: " + str(end - start))
    
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

    start = time.time()
    run_shell_command(samtools_host_subtract_command)
    end = time.time()
    print("Extracting sequences from the sam file from the snap step (after sorting it) took: " + str(end - start))

    #################################### SORTMERNA ############################
    aligned = dirpath + "/aligned"
    noRna = dirpath + "/noRrna"
    noRna1 = dirpath + "/noRrna_fwd.fq"
    noRna2 = dirpath + "/noRrna_rev.fq"

    # sortmerna index files are generated in the home directory by default, 
    # need to delete this as sortmerna will fail if this directory exists
    if os.path.exists(os.path.expanduser('~') + "/sortmerna"):
        shutil.rmtree(os.path.expanduser('~') + "/sortmerna")

    # avoid overlaps
    if os.path.isdir(dirpath + "/sortmerna"):
        shutil.rmtree(dirpath + "/sortmerna")

    sortmerna_command = "sortmerna" +\
                    " --ref " + args.sortmerna_rrna_database +\
                    " --aligned " + aligned +\
                    " --other " + noRna +\
                    " --fastx " +\
                    " --reads " + host_subtract_1 + " --reads " + host_subtract_2 +\
                    " --threads " + str(args.threads) +\
		            " --out2 TRUE " +\
		            " --paired_in TRUE" +\
                    " --workdir " + dirpath + "/sortmerna" # default path is home, if more than one workdir folder exists, the program will fail

    start = time.time()
    run_shell_command(sortmerna_command)
    end = time.time()
    print("Sortmerna rRNA depletion took: " + str(end - start))

    #################################### CLUMPIFY DEDUP #######################

    fullyQc1 = dirpath + "/fullyQc_fwd.fq"
    fullyQc2 = dirpath + "/fullyQc_rev.fq"

    clumpify_command = "clumpify.sh  in1=" + noRna1 + " in2=" + noRna2 +\
                            " out1=" + fullyQc1 + " out2=" + fullyQc2 + " dedupe=t &> " + dirpath + "/clumpify.log" + " -Xmx" + str(args.memory) + "g"
    
    start = time.time()
    run_shell_command(clumpify_command)
    end = time.time()
    print("Deduplication via clumpify took: " + str(end - start))

    # QUICK ALIGNMENT, JUST ALIGN REMAINING READS USING KRAKEN AGAINST KRAKEN_PLUS
    '''
    num_reads_bytes = subprocess.run(['grep', '-c', '.*', fullyQc1], stdout=subprocess.PIPE)
    num_reads_str = num_reads_bytes.stdout.decode('utf-8')
    # num_reads = int(num_reads_str.replace('', ''))/4 finally in int format, dividing by 4 because its in fastq format
    fast_mode_output = dirpath + "/fast_mode_output"
    kraken_command = "kraken2 --db " + args.kraken_db + " --threads " + str(args.threads) +\
                        " --output " + fast_mode_output + " --paired " + fullyQc1 + " " + fullyQc2
    
    run_shell_command(kraken_command)    
    kraken_res_out = dirpath + "/kraken_res_out"
    process_fast_mode_output(fast_mode_output, kraken_res_out, num_reads)
    fastAbundances = dirpath + "/fastAbundances.txt"
    get_lineage_info(kraken_res_out, fastAbundances, args.taxdump_location)
    fastAbundancesKrona = dirpath + "/fastAbundancesKrona.html"
    run_shell_command("ImportText.pl " + fastAbundances + " -o " + fastAbundancesKrona)
    '''
    
    #################### MEGAHIT ###########################
    megahit_out_path = args.dirpath

    if megahit_out_path[-1] == "/":
        megahit_out_path = megahit_out_path[0:-1]
    megahit_out_path = megahit_out_path + "/megahit_out"

    # deleting megahit output in the same location if it exists
    if os.path.isdir(megahit_out_path):
        shutil.rmtree(megahit_out_path)

    megahit_command = "megahit" + " -1 " + fullyQc1 +\
                        " -2 " + fullyQc2 +\
                        " -o " + megahit_out_path + " -t " + str(args.threads) # is an output directory 
    contigs = megahit_out_path + "/final_contigs.fa"

    start = time.time()
    run_shell_command(megahit_command)
    end = time.time()
    print("Assembly via megahit took: " + str(end - start))

    # need to check a few things to see if megahit failed by looking at the last line of the log file
    did_megahit_fail = False
    log_path = megahit_out_path + "/log"
    megahit_exit_message = subprocess.check_output('tail -n 1 ' + log_path, shell=True).decode('utf-8').strip('\n').replace("[","").replace("]","")

    if "Exit code" in megahit_exit_message:
        megahit_exit_status = megahit_exit_message.split(" ")[2]
        if megahit_exit_status.isnumeric() and int(megahit_exit_status) != 0:
            did_megahit_fail = True

    run_shell_command("mv " + megahit_out_path + "/final.contigs.fa " + contigs)

    if did_megahit_fail is False:
        # we must retrieve the unaligned reads
        reads_mapped_to_contigs_file_unsorted = dirpath + "/reads_mapped_to_contigs_unsorted.sam"
        align_reads_to_contigs_cmd = "bbwrap.sh" + " ref=" + contigs +\
                                    " in=" + fullyQc1 +\
                                    " in2=" + fullyQc2 +\
                                    " out=" + reads_mapped_to_contigs_file_unsorted + " -Xmx" + str(args.memory) + "g"
        
        start = time.time()
        run_shell_command(align_reads_to_contigs_cmd)
        end = time.time()
        print("bbwrap.sh alignment of the fullyQC reads to their contigs took: " + str(end - start))

        reads_mapped_to_contigs_file = dirpath + "/reads_mapped_to_contigs.sam"
        samtools_sort_command = "samtools sort -@ " + str(args.threads) + " -n " + reads_mapped_to_contigs_file_unsorted + " -o " + reads_mapped_to_contigs_file

        start = time.time()
        run_shell_command(samtools_sort_command)
        end = time.time()
        print("The time taken to sort the samtools file generated after mapping the reads onto the contigs is: " + str(end - start))

        # now lets retrieve the reads that did not align
        unassembled_reads_fwd = dirpath + "/unassembled_reads_fwd.fq"
        unassembled_reads_rev = dirpath + "/unassembled_reads_rev.fq"

        # same principle here as the host mapping step
        align_command = "samtools fastq -f 12 -1 " + unassembled_reads_fwd +\
                        " -2 " + unassembled_reads_rev + " " + reads_mapped_to_contigs_file
        
        start = time.time()
        run_shell_command(align_command)
        end = time.time()
        print("Retreiving the unassembled reads from the previous step took: " + str(end - start))

        # need to convert file above from fa to fq, simply done using seqtk
        contigs_fq = megahit_out_path + "/final_contigs.fq"
        seqtk_command = "seqtk" + " seq -F '#' " + contigs + " > " + contigs_fq
        run_shell_command(seqtk_command)

    else:
        unassembled_reads_fwd = dirpath + "/unassembled_reads_fwd.fq"
        unassembled_reads_rev = dirpath + "/unassembled_reads_rev.fq"
        shutil.copyfile(fullyQc1, unassembled_reads_fwd)
        shutil.copyfile(fullyQc2, unassembled_reads_rev)

    unassembled_reads_shorter_1 = dirpath + "/unassembled_reads_shorter_fwd.fq"
    unassembled_reads_shorter_2 = dirpath + "/unassembled_reads_shorter_rev.fq"
    unassembled_reads_longer_1 = dirpath + "/unassembled_reads_longer_fwd.fq"
    unassembled_reads_longer_2 = dirpath + "/unassembled_reads_longer_rev.fq"

    separate_reads_by_size(unassembled_reads_fwd, unassembled_reads_longer_1, unassembled_reads_shorter_1)
    separate_reads_by_size(unassembled_reads_rev, unassembled_reads_longer_2, unassembled_reads_shorter_2)

    # merge short reads, needed for blast alignment
    combined_file_sr_fq = dirpath + "/combined_sr_file.fq"
    merge_command = "seqtk mergepe " + unassembled_reads_shorter_1 + " " + unassembled_reads_shorter_2 + " > " + combined_file_sr_fq
    run_shell_command(merge_command)

    # convert above to fasta
    combined_file_sr_fa = dirpath + "/combined_sr_file.fa"
    convert_command = "seqtk seq -a " + combined_file_sr_fq + " > " + combined_file_sr_fa
    run_shell_command(convert_command)

    # need to merge paired end reads
    merged_pe = dirpath + "/merged_reads.fq"
    merge_command = "seqtk mergepe " + unassembled_reads_fwd + " " + unassembled_reads_rev + " > " + merged_pe
    run_shell_command(merge_command)
    
    combined_file_fq = dirpath + "/combined_reads_contigs_file.fq"

    if did_megahit_fail is False:
        # merge pe reads with contigs
        run_shell_command("cat " + contigs_fq + " > " + combined_file_fq)
        run_shell_command("cat " + merged_pe + " >> " + combined_file_fq)
    else:
        # if megahit failed, just runs merged_pe.
        run_shell_command("cat " + merged_pe + " > " + combined_file_fq)
        
    # need to convert above to fasta
    combined_file_fa = dirpath + "/combined_reads_contigs_file.fa"
    seqtk_command = "seqtk" + " seq -a " + combined_file_fq + " > " + combined_file_fa
    run_shell_command(seqtk_command)


    # creating alignment folder
    alignments_path = args.dirpath

    if alignments_path[-1] == "/":
        alignments_path = alignments_path[0:-1]
    
    alignments_path = alignments_path + "/alignments"

    if os.path.exists(alignments_path) is False:
        os.mkdir(alignments_path)

    # need to now go through each important output file to get relevant statistics for the output html

    summaryPath = dirpath + "/../summary"
    if os.path.exists(summaryPath) is False:
        os.mkdir(summaryPath)

    summaryFile = summaryPath + "/summary.txt"
    summaryFileWriter = open(summaryFile, "w")

    numReadsAtStart = 0
    numReadsAfterFastq = 0

    start = time.time()

    with open(fastp_json, "r") as fp_json:
        fastp_data = json.load(fp_json)
        numReadsAtStart = int(fastp_data["summary"]["before_filtering"]["total_reads"]/2)
        numReadsAfterFastq = int(fastp_data["summary"]["after_filtering"]["total_reads"]/2)

    summaryFileWriter.write("Start\t" + str(numReadsAtStart) + "\n")
    summaryFileWriter.write("Fastq\t" + str(numReadsAfterFastq) + "\n")
    summaryFileWriter.write("lowQuality\t" + str(numReadsAtStart - numReadsAfterFastq) + "\n")
    
    numReadsAfterStar = countNumSeqs(star_host_dedup1, False)
    summaryFileWriter.write("Star\t" + str(numReadsAfterStar) + "\n")

    if erccReadCounts is None:
        summaryFileWriter.write("erccReadCounts\t" + str("N/A") + "\n")
    else:
        summaryFileWriter.write("erccReadCounts\t" + str(erccReadCounts) + "\n")
    
    numReadsAfterSnap = countNumSeqs(host_subtract_1, False)
    summaryFileWriter.write("Snap\t" + str(numReadsAfterSnap) + "\n")

    if erccReadCounts is None:
        erccReadCounts = 0
    
    nonERCCHostReads = (numReadsAfterFastq - numReadsAfterStar - erccReadCounts) + (numReadsAfterStar - numReadsAfterSnap)

    summaryFileWriter.write("nonERCCHostReads\t" + str(nonERCCHostReads) + "\n")

    # go through clumpify log file
    with open(dirpath + "/clumpify.log", "r") as f:
        numReadsAfterSortmerna = 0
        nonHostRRNA = 0
        numDuplicates = 0
        for line in f:
            if "Reads In:" in line:
                numReadsAfterSortmerna = int(int(line.split(":")[1].strip())/2)
                nonHostRRNA = numReadsAfterSnap - numReadsAfterSortmerna
                summaryFileWriter.write("Sortmerna\t" + str(numReadsAfterSortmerna) + "\n")
                summaryFileWriter.write("nonHostRRNA\t" + str(nonHostRRNA) + "\n")
            if "Reads Out:" in line:
                numReadsAfterClumpify = int(int(line.split(":")[1].strip())/2)
                numDuplicates = numReadsAfterSortmerna - numReadsAfterClumpify
                summaryFileWriter.write("Clumpify\t" + str(numReadsAfterClumpify) + "\n")
                summaryFileWriter.write("duplicates\t" + str(numDuplicates) + "\n")

    totalReads = (numReadsAtStart - numReadsAfterFastq) + nonERCCHostReads + erccReadCounts + nonHostRRNA + numDuplicates + numReadsAfterClumpify
    summaryFileWriter.write("totalReadsChecked\t" + str(totalReads) + "\n")

    if did_megahit_fail is False:
        assemblyStats = subprocess.check_output('tail -n 2 ' + log_path + ' | head -n 1', shell=True).decode('utf-8').strip('\n').split(' ')

        numContigs = assemblyStats[2]
        summaryFileWriter.write("numContigs\t" + str(numContigs) + "\n")

        totalBases = assemblyStats[5]
        summaryFileWriter.write("totalContigsBases\t" + str(totalBases) + "\n")

        shortestContig = assemblyStats[8]
        summaryFileWriter.write("shortestContig\t" + str(shortestContig) + "\n")
        
        longestContig = assemblyStats[11]
        summaryFileWriter.write("longestContig\t" + str(longestContig) + "\n")

        avgContigLength = assemblyStats[14]
        summaryFileWriter.write("avgContigLength\t" + str(avgContigLength) + "\n")

        n50 = assemblyStats[17]
        summaryFileWriter.write("n50\t" + str(n50) + "\n")
    else:
        summaryFileWriter.write("numContigs\t" + "N/A" + "\n")
        summaryFileWriter.write("totalContigsBases\t" + "N/A" + "\n")
        summaryFileWriter.write("shortestContig\t" + "N/A" + "\n")
        summaryFileWriter.write("longestContig\t" + "N/A" + "\n")
        summaryFileWriter.write("avgContigLength\t" + "N/A" + "\n")
        summaryFileWriter.write("n50\t" + "N/A" + "\n")

    # will retreive numAssembledContigs at the finalisation step
    numLongReadsUnassembled = countNumSeqs(unassembled_reads_longer_1, False)
    numShortReadsUnassembled = countNumSeqs(unassembled_reads_shorter_1, False)

    summaryFileWriter.write("numReadsUnassembled\t" + str(int(numLongReadsUnassembled) + int(numShortReadsUnassembled)) + "\n")
    summaryFileWriter.write("numLongReadsUnassembled\t" + str(numLongReadsUnassembled) + "\n")
    summaryFileWriter.write("numShortReadsUnassembled\t" + str(numShortReadsUnassembled) + "\n")

    summaryFileWriter.close()

    end = time.time()
    print("Writing the statistics for the pipeline summary took: " + str(end - start))

    # zipping any unecessary files with pigz (supports multithreaded zipping)
    # no need to check for errors
    start = time.time()
    for toZip in [qc1, qc2, star_host_dedup1, star_host_dedup2, host_subtract_1, host_subtract_2, host_spare, 
                  snap_host_mapping, noRna1, noRna2, dirpath + "/star_host_ReadsPerGene.out.tab", dirpath + "/star_host_SJ.out.tab", 
                  aligned + "_fwd.fq", aligned + "_rev.fq"]:
        if os.path.exists(toZip + ".gz"):
            os.remove(toZip + ".gz")
        zip_command = "pigz " + toZip + " -p " + str(args.threads)
        run_shell_command(zip_command)
    
    end = time.time()
    print("Zipping the output files took: " + str(end - start))

    try:
        os.remove(snap_host_mapping_unsorted)
        os.remove(reads_mapped_to_contigs_file_unsorted)
    except:
        pass