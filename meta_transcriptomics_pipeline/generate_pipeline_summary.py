htmlCode = """
<html>
    <head>
        <title>Pipeline Steps</title>
        <style>
            body {
                display: flex;
                font-family: Helvetica;
                margin: 0;
                flex-direction: column;
                min-height: 100%; /* to ensure we fill the remaining space*/
            }

            .title {
                border: 5px;
                width: 100%;
                border-bottom: 1px gray;
            }
            
            .title_text {
                position: relative;
                left:15px;
            }

            .code_style {
                background: grey;
                box-sizing: border-box;
                padding: 10px;
                filter: brightness(182%);
            }

            .main_body {
                background: grey;
                filter: brightness(182%);
                width: 100%;
                height: 50%;
                overflow: auto;
                white-space: nowrap;
                display:flex;
                flex-direction: row;
                align-items: center;
                justify-content: space-between;
                flex-grow: 1; /* used to fill the remaining space in the body */
            }

            .box {
                background: rgba(64, 224, 208, 1.0);
                border-radius: 25px; /* for the curved borders */
                filter: brightness(50%);
                width: 200px;
                height: 100px;
                margin: 50px;
                flex-shrink: 0;
                display: flex;
                justify-content: center;
                align-items: center; /* although it doesnt appear so flex is good for positioning stuff in the center */
            }

            .bigger_box {
                background: rgba(252, 155, 171, 0.15);
                border-radius: 25px; /* for the curved borders */
                width: 300px;
                height: 400px;
                margin: 50px;
                flex-shrink: 0;
                display: flex;
                flex-direction: column;
                justify-content: center;
                align-items: center; /* although it doesnt appear so flex is good for positioning stuff in the center */
            }

            .box_text {
                color: white;
                position: relative;
                right: 0px;
                font-size: 25;
                text-align: center; /* so the text does not overflow to the right */
                white-space: initial;
                cursor: default; /* keep mouse cursor normal */
            }

            .box:hover {
                filter: brightness(40%);
            }

            .box:not(#nt_alignment):not(#nr_alignment):not(#get_read_abundances):not(#get_tpm_abundances):not(:last-child):after {
                content: "";
                height: 1px;
                background: black;
                width: 90px;
                position: absolute;
                right: -90px;
                top: 50%;
            }

            .box:not(#nt_alignment):not(#nr_alignment):not(#get_read_abundances):not(#get_tpm_abundances):not(:last-child):before {
                content: "";
                position: absolute;
                width: 0;
                height: 0;
                top: 50%;
                border-style: solid;
                border-width: 7px 0 7px 20px;
                border-color: transparent transparent transparent black;
                right: -100px;
                transform: translateY(-50%);
            }

            .bigger_box_line {
                content: "";
                height: 1px;
                background: black;
                width: 90px;
                position: absolute;
                left: -100px;
                top: 50%;
            }

            .bigger_box_arrow {
                content: "";
                position: absolute;
                width: 0;
                height: 0;
                top: 50%;
                border-style: solid;
                border-width: 7px 0 7px 20px;
                border-color: transparent transparent transparent black;
                transform: translateY(-50%);
                left: -20px;
                z-index: 10;
            }

            .popup_container {
                background: white;
                position: absolute;
                width: 50%;
                height: 100vh;
                left: -50%;
                overflow: auto;
                transition: transform 1s;
            }

            .popup_text {
                position: absolute;
                background: white;
                width: 100%;
                height: 100vh;
                text-align: left;
                box-sizing: border-box;
                padding: 20px;
            }

            .info {
                display: none;
            }
        </style>
    </head>
    <body>
        <div class="title">
            <div class="title_text">
                <h1>Pipeline Steps</h1>
            </div>
        </div>

        <div class="main_body">
            <div class="box" id="start">
                <div class="box_text">
                    Start
                </div>
            </div>
            <div class="box" id="qc">
                <div class="box_text">
                    QC
                </div>
            </div>
            <div class="box" id="human_subtraction">
                <div class="box_text">
                    Human Subtraction
                </div>
            </div>
            <div class="box" id="rrna_depletion">
                <div class="box_text">
                    rRNA Depletion
                </div>
            </div>
            <div class="box" id="deduplication">
                <div class="box_text">
                    Deduplication
                </div>
            </div>
            <div class="box" id="assembly">
                <div class="box_text">
                    Assembly
                </div>
            </div>
            <div class="bigger_box" id="alignments">
                <div class="box" id="nt_alignment">
                    <div class="box_text">
                        NT Alignment
                    </div>
                </div>
                <div class="box" id="nr_alignment">
                    <div class="box_text">
                        NR Alignment
                    </div>
                </div>
            </div>
            <div class="box" id="identify_best_alignments">
                <div class="bigger_box_line"></div>
                <div class="bigger_box_arrow"></div>
                <div class="box_text">
                    Identify Best Hits
                </div>
            </div>
            <div class="box" id="obtain_taxonomies">
                <div class="box_text">
                    Obtain TaxIds
                </div>
            </div>
            <div class="box" id="get_mean_alignment_scores">
                <div class="box_text">
                    Get Mean Alignment Scores
                </div>
            </div>
            <div class="box" id="get_taxid_lineages">
                <div class="box_text">
                    Get TaxId Lineages
                </div>
            </div>
            <div class="bigger_box" id="abundances">
                <div class="box" id="get_read_abundances">
                    <div class="box_text">
                        Get Read Abundances
                    </div>
                </div>
                <div class="box" id="get_tpm_abundances">
                    <div class="box_text">
                        Get TPM Abundances
                    </div>
                </div>
            </div>
        </div>
        <div class="popup_container">
            <div class="popup_text">
                Lorem Ipsum
            </div>
        </div>
        <script>
            var popup_container = document.getElementsByClassName("popup_container")[0]
            var popup = document.getElementsByClassName("popup_text")[0]
            var boxes = document.getElementsByClassName("box")

            function clicked({ target }) {
                elem_id = target.id;
                if (target.className === "box_text") {
                    elem_id = target.parentNode.id
                }
                popup_container.style.transform = "translateX(100%)";
                var elem = document.getElementById(elem_id + "_info");
                popup.innerHTML = elem.innerHTML;
            }

            /* from: https://htmldom.dev/check-if-an-element-is-a-descendant-of-another/ */

            const isDescendant = function (parent, child) {
                let node = child.parentNode;
                while (node) {
                    if (node === parent) {
                        return true;
                    }
                    // Traverse up to the parent
                    node = node.parentNode;
                }
                // Go up until the root but couldn't find the `parent`
                return false;
            };

            window.addEventListener("click", ({ target }) => {
                if (target.matches(".box") || target.matches(".box_text")) return; // checking item classname

                if (target.isEqualNode(popup_container) === true || isDescendant(popup_container, target) === true) return;
                popup_container.style.transform = "translateX(-100%)";
            });

            for (let i = 0; i < boxes.length; i++) {
                boxes[i].addEventListener("click", clicked)
            }
        </script>
    </body>
    <div class="info" id="start_info">
        <h1 style="text-align: center;">Start</h1> <br />
        <b>Num Reads at Start:</b> <br />
    </div>
    <div class="info" id="qc_info">
        <h1 style="text-align: center;">Quality Control</h1>
        <p>
            <b>Description:</b> This step uses Fastp to remove reads that have low quality, low complexity or are too short (less than 50bp). It also strips away the adapters.
        </p>
        <b>Num Reads After Fastq:</b> <br>
        <b>Outputs:</b> fastp_1.fastq, fastp_2.fastq <br> <br>
        <div class="code_style">
            <code>
                fastp --in1 input1.fastq.gz --in2 input2.fastq.gz \​ <br>
                --out1 fastp_1.fastq --out2 fastp_1.fastq \​ <br>
                -b 100 -B 100 \​ <br>
                --qualified_quality_phred 15 \​ <br>
                --unqualified_percent_limit 45 \​ <br>
                --length_required 50 \​ <br>
                --low_complexity_filter \​ <br>
                --detect_adapter_for_pe \​ <br>
                --thread N <br>
            </code>
        </div>
        <p>
            <b>Explanation of Code:</b> -b and -B will ensure both pairs will not have more than 100 bases by trimming away any excess bases. --qualified_qualifty_phred is the quality value (Q) that a base must have for it to be "qualified". 15 means we want a quality score of 15 (i.e. Q15). If the percentage of unqualified bases in a read exceeds --unqualified_percent_limit, the read will be discarded.
        </p>
    </div>
    <div class="info" id="human_subtraction_info">
        <h1 style="text-align: center;">Human Subtraction</h1>
        <b>Description:</b> This step removes human sequences from the library. It uses STAR to align the library against an index of the HG38 assembly built using STAR and we extract the reads that failed to map (i.e. the non-human reads). To maximise the amount of reads removed, this step is repeated using SNAP to align the remaining reads against the HG38 index built using SNAP. <br> <br>
        <b>Num Reads Remaining After Human Depletion:</b> <br>
        <b>Num Reads Remaining After STAR Command:</b> <br>
        <b>Num Reads Remaining After SNAP Command:</b> <br>
        <b>Inputs:</b> fastp_1.fastq, fastp_2.fastq <br>
        <b>Outputs:</b> host_depleted1.fastq, host_depleted2.fastq <br> <br>
        <div class="code_style">
            <code>
                # STAR ALIGNMENT <br>
                # The STAR step outputs the reads that failed to map (the nonhuman reads) into star_host_Unmapped_1/2.fastq<br> <br>
                STAR --genomeDir (star_host_index) \​ <br>
                --readFilesIn fastp_1.fastq fastp_2.fastq \​ <br>
                --outFileNamePrefix "/star_host_" \​ <br>
                --outFilterMultimapNmax 99999 \​ <br>
                --outFilterScoreMinOverLread 0.5 \​ <br>
                --outFilterMatchNminOverLread 0.5 \​ <br>
                --outFilterMismatchNmax 999 \​ <br>
                --outSAMtype BAM SortedByCoordinate \​ <br>
                --outReadsUnmapped Fastx" \​ <br>
                --outSAMattributes Standard \​ <br>
                --quantMode TranscriptomeSAM GeneCounts \​ <br>
                --clip3pNbases 0 \​ <br>
                --runThreadN N <br> <br>

                # SNAP ALIGNMENT <br>
                snap-aligner paired (snap_host_index) \ ​<br>
                star_host_Unmapped_1.fastq star_host_Unmapped_2.fastq \​ <br>
                 -o snap_host_mapped.bam -t N + -I <br> <br>

                # USE SAMTOOLS TO RETREIVE UNALIGNED READS FROM SNAP OUTPUT <br> <br>
                samtools fastq -f 12 -@ N \​ <br>
                -1 host_depleted1.fastq \​ <br>
                -2 host_depleted2.fastq \​ <br>
                -s undetermined.fastq \​ <br>
                snap_host_mapped.bam <br> <br>
            </code>
        </div>
        <p>
            <b>Explanation of Code:</b> In samtools fastq, -f 12 means we want to use filtering option 12 to remove reads where both pairs failed to map. This means we want to filter out reads that failed to align (4) and mates that failed to map (8). Hence 4 + 8 = 12.
        </p>
    </div>
    <div class="info" id="rrna_depletion_info">
        <h1 style="text-align: center;">rRNA Depletion</h1>
        <b>Description:</b> This step removed any rRNA reads from the library by aligning them against the "default" sortmerna database available from <a href="https://github.com/sortmerna/sortmerna/releases">here</a>. It is important to run this step after the human read depletion step as sortmerna is very slow and RAM hungry, hence removing as much reads as possible prior to this step will increase performance. <br> <br>
        <b>Num Reads After rRNA Depletion:</b> <br>
        <b>Inputs:</b> host_depleted1.fastq, host_depleted2.fastq  <br>
        <b>Outputs:</b> noRna_fwd.fq, noRna_rev.fq <br> <br>
        <div class="code_style">
            <code>
                sortmerna --ref (sortmerna_rrna_database) \​ <br>
                --fastx \​ <br>
                --aligned aligned \​ <br>
                --other noRrna \​ <br>
                --reads host_depleted1.fastq --reads host_depleted2.fastq \​ <br>             
                --out2 TRUE \​ <br>
                --paired_in TRUE \​ <br>
                --threads N \ <br>
            </code>
        </div>
    </div>
    <div class="info" id="deduplication_info">
        <h1 style="text-align: center;">Deduplication</h1>
        <p>
            <b>Description:</b> Remove duplicate sequences, these may have been the same read amplified multiple times during PCR. This is the final step in the quality control phase of the pipeline and the reads we get out from this step will be used for quantification and analysis.
        </p>
        <b>Num Reads After Deduplication:</b> <br>
        <b>Inputs:</b> noRna_fwd.fq, noRna_rev.fq  <br>
        <b>Outputs:</b> fullyQc_fwd.fq, fullyQc_rev.fq <br> <br>
        <div class="code_style">
            <code>
                clumpify.sh in1=noRna_fwd.fq in2=noRna_rev.fq \​ <br>
                out1=fullyQc_fwd.fq out2=fullyQc_rev.fq dedupe=t​ <br>
            </code>
        </div>
    </div>
    <div class="info" id="assembly_info">
        <h1 style="text-align: center;">Assembly</h1>
        <p>
            <b>Description:</b> Assemble our reads using megahit and identify unassembled/assembled reads.
        </p>
        <b>Num Contigs Assembled:</b> <br>
        Num Reads Assembled into Contigs:</b> <br>
        <b>Num Reads Unassembled:</b> <br>
        <b>Total Contig Bases:</b> <br>
        <b>Shortest Contig:</b> <br>
        <b>Longest Contig:</b> <br>
        <b>Average Contig Length:</b> <br>
        <b>N50 Score:</b> <br>
        <b>Number of Unassembled Short Reads:</b> <br>
        <b>Number of Unassembled Long Reads:</b> <br>
        <b>Inputs:</b> fullyQc_fwd.fq, fullyQc_rev.fq <br>
        <b>Outputs:</b> megahit_out/final_contigs.fa, unassembled_reads_fwd.fq, unassembled_reads_rev.fq, reads_belonging_to_contigs.txt <br> <br>
        <div class="code_style">
            <code>
                megahit -1 fullyQc_fwd.fq -2 fullyQc_rev.fq \​ <br>
                -o + megahit_out -t N <br> <br> 

                # map reads onto the assembly again (megahit does not provide this)​ <br>
                bbwrap.sh ref=megahit_out/final_contigs.fa \​ <br>
                in=fullyQc_fwd.fq in2=fullyQc_rev.fq \​ <br>
                -out=reads_mapped_to_contigs.sam <br> <br>

                # use identical command to the human mapping step to retreive unaligned reads <br>
                samtools fastq -f 12 \​ <br>
                -1 unassembled_reads_fwd.fq -2 unassembled_reads_rev.fq \​ <br>
                reads_mapped_to_contigs.sam <br> <br>

                # aside from the above, we also go through reads_mapped_to_contigs.sam once again to extract the assembled reads. <br>
                # for each query, we select the hit with the lowest edit distance but if multiple hits share this value, we select the hit among these with the highest bitscore. <br>
                # the reads_belonging_to_contigs.txt file contains the read name and its corresponding contig, sorted by contig name. <br>
            </code>
        </div>
    </div>
    <div class="info" id="nt_alignment_info">
        <h1 style="text-align: center;">NT BLAST Alignment</h1>
        <p>
            <b>Description:</b> Aligns the unassembled reads and contigs against the NCBI NT database using BLAST.
        </p>
        <b>Num Unique BLAST Hits:</b> <br>
        <b>Inputs:</b> unassembled_reads_fwd.fq, unassembled_reads_rev.fq, megahit_out/final_contigs.fq (as combined_file.fa where the paired end reads are interleaved)  <br>
        <b>Outputs:</b> nt_alignments_file.tsv <br> <br>
        <div class="code_style">
            <code>
                blastn –task megablast \ <br> 
                -query combined_file.fa \ <br> 
                -db nt \ <br> 
                -out nt_alignments_file.tsv \​ <br>
                -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \ <br>
                -max_target_seqs 10 \ <br>
                -num_threads N <br>
            </code>
        </div>
    </div>
    <div class="info" id="nr_alignment_info">
        <h1 style="text-align: center;">NR BLAST Alignment</h1>
        <p>
            <b>Description:</b> Use diamond to align the unassembled reads and contigs against an index of the NCBI NR database built using diamond.
        </p>
        <b>Num Unique Diamond Hits:</b> <br>
        <b>Inputs:</b> unassembled_reads_fwd.fq, unassembled_reads_rev.fq, megahit_out/final_contigs.fq (as combined_file.fq where the paired end reads are interleaved)  <br>
        <b>Outputs:</b> nr_alignments_file.tsv <br> <br>
        <div class="code_style">
            <code>
                diamond blastx \ <br>
                --db nr_index \ <br>
                --query combined_file.fq \ <br>
                --out nr_alignments_file.tsv \​ <br>
                --mid-sensitive \ <br>
                --max-target-seqs 10 \ <br> 
                --masking 0 –c 1 –b 6 \​ <br>
                --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen \ <br>
                --threads N <br>
            </code>
        </div>
        <p>
            <b>Explanation of Code:</b> --mid-sensitive refers to the sensitivity option, --masking 0 disables repeat masking of the query and reference sequences, -b refers to the number of sequence letters to load in memory at any one time (in billions) and -c refers to the number of chunks to divide each seed index when processing it.
        </p>
    </div>
    <div class="info" id="identify_best_alignments_info">
        <h1 style="text-align: center;">Identify Best Hits</h1>
        <p>
            <b>Description:</b> Selects the best hits from the alignment stage. 
        </p>
        <b>Number of Best Hits from the NT Alignment:</b> <br>
        <b>Number of Best Hits from the NR Alignment:</b> <br>
        <b>Number of Unassembled Reads Successfully Aligned:</b> <br>
        <b>Number of Unassembled Contigs Successfully Aligned:</b> <br>
        <b>Number of Unique Accessions</b> <br>
        <b>Inputs:</b> nt_alignments_file.tsv, nr_alignments_file.tsv  <br>
        <b>Outputs:</b> best_nt_scores.tsv, best_nr_scores.tsv <br>
        <p>
            <b>Explanation of Step:</b> Combine the nt and nr alignments from the previous step to identify the best hits for each sequence. Hits are selected based on their e-value and if multiple hits have the same e-value, we select the hit from these with the highest bitscore. Where the hit will be placed will depend on whether the hit originated from the NT or NR alignment. The final files will store the read/contig name, its length (in the final column), the accession of the sequence it mapped onto, and the alignment's e-value, bitscore and percentage identity.
        </p>
    </div>
    <div class="info" id="obtain_taxonomies_info">
        <h1 style="text-align: center;">Obtain Taxonomic IDs</h1>
        <p>
            <b>Description:</b> Obtain the taxids of the accessions, and hence the taxids of the contigs and the assembled/unassembled reads. 
        </p>
        <b>Number of Mapped Accessions:</b> <br>
        <b>Number of Unique Taxids:</b> <br>
        <b>Number of Contigs with a Taxid:</b> <br>
        <b>Number of Unassembled Reads with a Taxid</b> <br>
        <b>Number of Assembled Reads with a Taxid</b> <br>
        <b>Outputs:</b> contigs_reads_accessions_taxids.txt, assembled_reads_taxids_temp.txt<br>
        <p>
            <b>Explanation of Step:</b> For both best_nt_scores.tsv and best_nr_scores.tsv, we sort them individually by their accessions and extract their query and accession IDs (ignoring version number). We then run them against one or more corresponding mapping files (which are sorted by accession ID) and this will take O(N) time (where N refers to the number of records in the mapping file). We will compare the accessions from best_nt_scores.tsv against dead_wgs.accession2taxid, nucl_gb.accession2taxid, nucl_wgs.accession2taxid and nucl_wgs.accession2taxid.EXTRA and we will compare the accessions from best_nr_scores.tsv against prot.accession2taxid.
            These files are available from <a href="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/">NCBI</a> and they are already sorted by their accessions. The file outputted in this step is contigs_reads_accessions_taxids.txt, which contains the query ID and its corresponding accession and taxonomic id, sorted by query ID.
            If you notice a difference in the number of unassembled reads/contigs with a taxid here when compared to the previous step, this is because the accession of the missing sequences are not present in any databases searched in this step. This section was written in Cython, as the accession to taxid mapping files are too large to process with Python. <br> <br>
            We also need to pair the assembled reads with the taxid of its parent contig, we can simply do this by fetching the taxid of the contigs in contigs_reads_accessions_taxids.txt and map them to the contigs in reads_belonging_to_contigs.tsv. The file assembled_reads_taxids_temp.txt will contain the read name and its corresponding taxid. 
        </p>
    </div>
    <div class="info" id="get_mean_alignment_scores_info">
        <h1 style="text-align: center;">Get Mean Alignment Scores</h1>
        <p>
            <b>Description:</b> Get the average alignment scores for each taxid.
        <p>
            <b>Explanation of Step:</b> Combine the results from best_nt_scores.tsv and best_nr_score.tsv and match their accessions to their taxid from contigs_reads_accessions_taxids.txt, which in turn is converted to its Scientific name.
            Final file is called species_avg_alignment_scores.txt, which contains the corresponding scientific names for its query, followed by the e-value, bitscore, percentage identity and query length of that hit. Finally, we combine all identical scientific names and get the mean of these scores and query length.
        </p>
    </div>
    <div class="info" id="get_taxid_lineages_info">
        <h1 style="text-align: center;">Get TaxID Lineages</h1>
        <p>
            <b>Description:</b> Obtain the lineages of the unique taxids from contigs_reads_accessions_taxids.txt and convert the taxids to their scientific names. <br>
            <b>Number of Taxid Lineages Successfully Identified:</b> <br>
            <b>Number of Taxid Lineages Not Identified:</b> <br>
        <p>
            <b>Explanation of Step:</b> For each taxid in contigs_reads_accessions_taxids.txt, we use the Python library ete3 to fetch the parent taxids of its lineages for these ranks [superkingdom, phylum, class, order, family and genus]. If for whatever reason any parent taxid in our lineage does not exist, call that rank "Unknown".
            Also, if our current taxid is higher than rank "species", discard this taxid. 
            We now convert all unique taxids to their scientific names, using the corresponding taxid to scientific name mappings in names.dmp from taxdump from NCBI
            We did the same thing to map our taxids to its Scientific names in the previous step.
        </p>
    </div>
    <div class="info" id="get_read_abundances_info">
        <h1 style="text-align: center;">Get Read Abundances</h1>
        <p>
            <b>Description:</b> Get the abundance of a sample using a simple read counts method of all reads with a taxid.
        </p>
        <b>Inputs:</b> contigs_reads_accessions_taxids.txt, assembled_reads_taxids_temp.txt <br>
        <b>Outputs:</b> abundanceReadMethod.txt <br>
        <p>
            <b>Explanation of Step:</b> Count the number of assembled/unassembled reads that map onto each taxid in contigs_reads_accessions_taxids.txt and assembled_reads_taxids_temp.txt and then map the taxids onto their lineages. The final file abundanceReadMethod.txt contains the number of reads mapped on a taxid (as a percentage of all reads after the deduplication step) and the lineage of the taxid.
        </p>
    </div>
    <div class="info" id="get_tpm_abundances_info">
        <h1 style="text-align: center;">Get TPM Abundances</h1>
        <p>
            <b>Description:</b> Get the abundance using a more complicated TPM method that normalises for "gene" length and sequencing depth.
        </p>
        <b>Inputs:</b> contigs_reads_accessions_taxids.txt, assembled_reads_taxids_temp.txt <br>
        <b>Outputs:</b> tpmAbundancesKrona.html <br>
        <p>
            <b>Explanation of Step:</b> Due to the lack of gene annotation data, this step uses the RSEM method that treats the contigs as genes. We also treat unassembled reads as "contigs", where they will be treated as having 1 read mapped onto them (i.e. the read maps onto itself) and it will have a length equivalent to the N50 score of the assembly (to prevent a calculation bias in favour of the shorter reads). <br> <br>
            We will count the number of reads assigned to each contig in reads_belonging_to_contigs.txt. We will also obtain the length of all contigs and "length" of all unassembled reads (they are assigned the N50 score of the assembly) but only if they have a taxid mapped onto them. For each "contig", we will place the "contig", their lengths as well as the number of reads mapped onto the contig and the taxid that they map onto in contig_unaligned_read_counts_len_taxid.txt.  <br> <br>
            From this file, we will calculate the TPM scores for each taxid and then map the taxids onto their lineage. The final file tpmAbundancesKrona.html  will contain the TPM score for each taxid and the lineages of the taxid.
        </p>
    </div>
</html>
"""

def generate_pipeline_summary(summaryFile, outputFile):
    global htmlCode
    with open(summaryFile, "r") as f:
        for record in f:
            curr = record.split("\t")
            if len(curr) < 2:
                continue
            curr[1] = curr[1].strip()
            if curr[0] == "Start":
                htmlCode = htmlCode.replace("Num Reads at Start:</b> ", "Num Reads at Start:</b> " + str(curr[1]))
            if curr[0] == "Fastq":
                htmlCode = htmlCode.replace("Num Reads After Fastq:</b> ", "Num Reads After Fastq:</b> " + str(curr[1]))
            if curr[0] == "Snap":
                htmlCode = htmlCode.replace("Num Reads Remaining After Human Depletion:</b> ", "Num Reads Remaining After Human Depletion:</b> " + str(curr[1]))
                htmlCode = htmlCode.replace("Num Reads Remaining After SNAP Command:</b> ", "Num Reads Remaining After SNAP Command:</b> " + str(curr[1]))
            if curr[0] == "Star":
                htmlCode = htmlCode.replace("Num Reads Remaining After STAR Command:</b> ", "Num Reads Remaining After STAR Command:</b> " + str(curr[1]))
            if curr[0] == "Sortmerna":
                htmlCode = htmlCode.replace("Num Reads After rRNA Depletion:</b> ", "Num Reads After rRNA Depletion:</b> " + str(curr[1]))
            if curr[0] == "Clumpify":
                htmlCode = htmlCode.replace("Num Reads After Deduplication:</b> ", "Num Reads After Deduplication:</b> " + str(curr[1]))
            if curr[0] == "numContigs":
                htmlCode = htmlCode.replace("Num Contigs Assembled:</b> ", "Num Contigs Assembled:</b> " + str(curr[1]))
            if curr[0] == "numAssembledReads":
                htmlCode = htmlCode.replace("Num Reads Assembled into Contigs:</b> ", "Num Reads Assembled into Contig:</b> " + str(curr[1]))
            if curr[0] == "numReadsUnassembled":
                htmlCode = htmlCode.replace("Num Reads Unassembled:</b> ", "Num Reads Unassembled:</b> " + str(curr[1]))
            if curr[0] == "totalContigsBases":
                htmlCode = htmlCode.replace("Total Contig Bases:</b> ", "Total Contig Bases:</b> " + str(curr[1]))
            if curr[0] == "shortestContig":
                htmlCode = htmlCode.replace("Shortest Contig:</b> ", "Shortest Contig:</b> " + str(curr[1]))
            if curr[0] == "longestContig":
                htmlCode = htmlCode.replace("Longest Contig:</b> ", "Longest Contig:</b> " + str(curr[1]))
            if curr[0] == "avgContigLength":
                htmlCode = htmlCode.replace("Average Contig Length:</b> ", "Average Contig Length:</b> " + str(curr[1]))
            if curr[0] == "n50":
                htmlCode = htmlCode.replace("N50 Score:</b> ", "N50 Score:</b> " + str(curr[1]))
            if curr[0] == "numShortReadsUnassembled":
                htmlCode = htmlCode.replace("Number of Unassembled Short Reads:</b> ", "Number of Unassembled Short Reads:</b> " + str(curr[1]))
            if curr[0] == "numLongReadsUnassembled":
                htmlCode = htmlCode.replace("Number of Unassembled Long Reads:</b> ", "Number of Unassembled Long Reads:</b> " + str(curr[1]))
            if curr[0] == "numUniqueNTHits":
                htmlCode = htmlCode.replace("Num Unique BLAST Hits:</b> ", "Num Unique BLAST Hits:</b> " + str(curr[1]))
            if curr[0] == "numUniqueNRHits":
                htmlCode = htmlCode.replace("Num Unique Diamond Hits:</b> ", "Num Unique Diamond Hits:</b> " + str(curr[1]))
            if curr[0] == "numBestNTHits":
                htmlCode = htmlCode.replace("Number of Best Hits from the NT Alignment:</b> ", "Number of Best Hits from the NT Alignment:</b> " + str(curr[1]))
            if curr[0] == "numBestNRHits":
                htmlCode = htmlCode.replace("Number of Best Hits from the NR Alignment:</b> ", "Number of Best Hits from the NR Alignment:</b> " + str(curr[1]))
            if curr[0] == "numMappedUnassembledReads":
                htmlCode = htmlCode.replace("Number of Unassembled Reads Successfully Aligned:</b> ", "Number of Unassembled Reads Successfully Aligned:</b> " + str(curr[1]))
            if curr[0] == "numMappedContigs":
                htmlCode = htmlCode.replace("Number of Unassembled Contigs Successfully Aligned:</b> ", "Number of Unassembled Contigs Successfully Aligned:</b> " + str(curr[1]))
            if curr[0] == "totalUniqueAccessions":
                htmlCode = htmlCode.replace("Number of Unique Accessions</b> ", "Number of Unique Accessions</b> " + str(curr[1]))
            if curr[0] == "numMappedAccessions":
                htmlCode = htmlCode.replace("Number of Mapped Accessions:</b> ", "Number of Mapped Accessions:</b> " + str(curr[1]))
            if curr[0] == "numUniqueTaxids":
                htmlCode = htmlCode.replace("Number of Unique Taxids:</b> ", "Number of Unique Taxids:</b> " + str(curr[1]))
            if curr[0] == "numContigsWithTaxids":
                htmlCode = htmlCode.replace("Number of Contigs with a Taxid:</b> ", "Number of Contigs with a Taxid:</b> " + str(curr[1]))
            if curr[0] == "numUnassembledReadsWithTaxids":
                htmlCode = htmlCode.replace("Number of Unassembled Reads with a Taxid</b> ", "Number of Unassembled Reads with a Taxid</b> " + str(curr[1]))
            if curr[0] == "numAssembledReadsWithTaxids":
                htmlCode = htmlCode.replace("Number of Assembled Reads with a Taxid</b> ", "Number of Assembled Reads with a Taxid</b> " + str(curr[1]))
            if curr[0] == "numTaxidsWithRankSpecies":
                htmlCode = htmlCode.replace("Number of Taxid Lineages Successfully Identified:</b> ", "Number of Taxid Lineages Successfully Identified:</b> " + str(curr[1]))
            if curr[0] == "numTaxidsWithoutRankSpecies":
                htmlCode = htmlCode.replace("Number of Taxid Lineages Not Identified:</b> ", "Number of Taxid Lineages Not Identified:</b> " + str(curr[1]))

    wf = open(outputFile, "w")
    wf.write(htmlCode)
    wf.close()