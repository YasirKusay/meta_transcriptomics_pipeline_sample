import argparse
import multiprocessing
from meta_transcriptomics_pipeline.blast_sr import blast_sr
from meta_transcriptomics_pipeline.diamond import diamond
from meta_transcriptomics_pipeline.finalisation import finalisation
from meta_transcriptomics_pipeline.minimap2_contigs import minimap2_contigs
from meta_transcriptomics_pipeline.minimap2_reads import minimap2_reads
from meta_transcriptomics_pipeline.preprocessing import preprocessing

def parse_args():
    parser = argparse.ArgumentParser(
        description="Detect contamination in an assembly"
    )

    subparsers = parser.add_subparsers()

    parser_preprocessing = subparsers.add_parser(
        "preprocessing", help="Perform the preprocessing step prior to the alignment on the input sequences including quality check, assembly and a preliminary prediction of the composition of the sample via kraken"
    )

    ################# PREPROCESSING ###################

    parser_preprocessing.add_argument(
        "inp1",
        type=str,
        help="Paired end read file 1 input path (must be in fastq format)"
    )

    parser_preprocessing.add_argument(
        "inp2",
        type=str,
        help="Paired end read file 2 input path (must be in fastq format)"
    )

    parser_preprocessing.add_argument(
        "sortmerna_rrna_database",
        type=str,
        help="Path to the sortmerna rRNA database used in the rRNA depletion step.\n" + \
        "Download the 'database.tar.gz' zipped folder from https://github.com/biocore/sortmerna/releases" + \
        "that matches your sortmerna version and select just a single database file from it.\n" + \
        "Our pipeline was tested using 'smr_v4.3_default_db.fasta'."
    )

    parser_preprocessing.add_argument(
        "star_host_index",
        type=str,
        help="Path to the pre-built STAR host index. The complete host genome needs to be indexed. For the human genome, we used GRCh38_latest_genomic.fna available from https://www.ncbi.nlm.nih.gov/genome/guide/human/ and index it via STAR.\n"
    )

    parser_preprocessing.add_argument(
        "snap_host_index",
        type=str,
        help="Path to the pre-built SNAP host index folder. The complete host genome needs to be indexed. For the human genome, we used GRCh38_latest_genomic.fna available from https://www.ncbi.nlm.nih.gov/genome/guide/human/ and index it via SNAP.\n" +\
        "SNAP generates a folder containing the index and the path to this folder needs to be specified."
    )

    parser_preprocessing.add_argument(
        "kraken_plus_db",
        type=str,
        help="The location of the kraken_plus database used to predict the composition of the sample prior to the alignment step. Available to download from https://benlangmead.github.io/aws-indexes/k2"
    )

    parser_preprocessing.add_argument(
        "dirpath",
        type=str,
        help="Path to a directory where to place the intermediate outputs. Needs to be the same for all subcommands"
    )

    parser_preprocessing.add_argument(
        "--qualified_quality_phred",
        type=str,
        default="15",
        help=" DURING QC. The quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified"
    )

    parser_preprocessing.add_argument(
        "--unqualified_percent_limit",
        type=str,
        default="45",
        help=" DURING QC. How many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%"
    )

    parser_preprocessing.add_argument(
        "--average_qual",
        type=str,
        default="0",
        help=" DURING QC. if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement (int [=0])"
    )

    parser_preprocessing.add_argument(
        "--length_required",
        type=str,
        default="50",
        help="DURING QC. Reads shorter than length_required will be discarded, default is 15. (int [=15])"
    )

    parser_preprocessing.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help="The number of cpu threads to use"
    )

    parser_preprocessing.set_defaults(func=preprocessing)

    parser_minimap_contigs = subparsers.add_parser(
        "minimap_contig_alignment", help="Perform the minimap alignment of the contigs generated in the previous step against the minimap2 index."
    )

    minimap2_index_message = "Path to the minimap2 index file to align reads/contigs. Needs to be built prior to this step using minimap2 index on the NCBI NT sequences (available from https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/). \n" + \
        "WARNING! PAY CAREFUL ATTENTION TO THE SIZE OF THE NCBI NT SEQUENCES AS IT IS CURRENTLY (AS OF FEBRUARY 2023) NEARLY 1TB AND THE MINIMAP2 INDEX WILL BE AT LEAST 2.5 TIMES THE SIZE OF THAT!"

    parser_minimap_contigs.add_argument(
        "minimap2_index",
        type=str,
        help=minimap2_index_message
    )

    parser_minimap_contigs.add_argument(
        "dirpath",
        type=str,
        help="Path to the directory used by the preprocessing step (which will store the contigs). Will place additional files here."
    )

    parser_minimap_contigs.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help="The number of cpu threads to use"
    )

    parser_minimap_contigs.set_defaults(func=minimap2_contigs)

    parser_minimap_reads = subparsers.add_parser(
        "minimap_long_read_alignment", help="Perform the minimap alignmnet against the minimap2 index using the reads (that are > 100bp) that failed to assemble in the pre-processing step."
    )

    parser_minimap_reads.add_argument(
        "minimap2_index",
        type=str,
        help=minimap2_index_message
    )

    dirpath_help = "Path to the directory used by the preprocessing step (which will store the unassembled reads). Will place additional files here."

    parser_minimap_reads.add_argument(
        "dirpath",
        type=str,
        help=dirpath_help
    )

    parser_minimap_reads.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help="The number of cpu threads to use"
    )

    parser_minimap_reads.set_defaults(func=minimap2_reads)

    parser_blast_reads = subparsers.add_parser(
        "blast_short_read_alignment", help="Perform the minimap alignmnet against the blast NCBI NT index using the reads (that are < 100bp) that failed to assemble in the pre-processing step. Unlike the former alignment, this uses the 'sr' (i.e. short read) option."
    )

    parser_blast_reads.add_argument(
        "blast_ncbi_nt_database",
        type=str,
        help="The path to the blast ncbi nt index. This needs to include not just the location, but the prepended name of the index files.\n" + \
        "Can simply be downloaded from https://ftp.ncbi.nlm.nih.gov/blast/db/ (fetching the files that have nt.[0-9]+"
    )

    parser_blast_reads.add_argument(
        "dirpath",
        type=str,
        help=dirpath_help
    )

    parser_blast_reads.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help="The number of cpu threads to use"
    )

    parser_blast_reads.set_defaults(func=blast_sr)

    parser_diamond = subparsers.add_parser(
        "diamond_alignment", help="Perform the diamond alignment on all contigs and unassembled reads against the diamond index built using NCBI NR."
    )

    parser_diamond.add_argument(
        "diamond_index",
        type=str,
        help="Path to the diamond index file to align reads/contigs. Needs to be built using nr.gz available from https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/ using diamond index."
    )

    parser_diamond.add_argument(
        "dirpath",
        type=str,
        help="Path to the directory used by the preprocessing step (which will store the unassembled reads and contigs). Will place additional files here."
    )

    parser_diamond.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help="the number of cpu threads to use"
    )

    parser_diamond.set_defaults(func=diamond)

    parser_finalisation = subparsers.add_parser(
        "finalisation", help="Perform taxnonomic identification, decontamination (optional) and quantification.\n" + \
        "Final outputs can be found in readAbundancesKrona.html which shows the Krona plot generated using the read count calculation and tpmAbundancesKrona.html which shows the Krona plot generated using the TPM calculation."
    )

    parser_finalisation.add_argument(
        "dirpath",
        type=str,
        help="Path to the directory used by the preprocessing step and alignment steps (which will store the alignment output). Places quantification data here."
    )

    parser_finalisation.add_argument(
        "taxdump_location",
        type=str,
        help="Path to the taxdump folder that contains information such as the lineages of the species. Available from: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/"
    )

    parser_finalisation.add_argument(
        "nucl_prot_accession_taxid_mapping_files_loc",
        type=str,
        help="Specify location that map sequence accessions (found in the NCBI NT/NR database) to their respective taxids. Files available from: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/. \n" + \
        "We used dead_nucl.accession2taxid.gz, dead_wgs.accession2taxid.gz , nucl_gb.accession2taxid.gz, nucl_wgs.accession2taxid.EXTRA.gz and nucl_wgs.accession2taxid.gz for nucleotide mappings and " + \
        " dead_prot.accession2taxid.gz and pdb.accession2taxid.gz for protein mappings. The file can automatically detect what kind of mapping file it is (either nucleotide accessions/protein accessions)."
    )

    parser_finalisation.add_argument(
        "--decontamination_input_files",
        type=str,
        help="The path to the place that stores kraken outputs of sequences from a similar environment as the input sequence. These files will be used in the decontamination step.\n" + \
            "This location has 2 folders: 'controls' where you will store the kraken output for the control sequences and 'others' where you will store the kraken output for the other sequences.\n" + \
            "If controls is empty or does not exist, we will be using the 'SQUEEGEE' program to identify contaminants (does not require controls), otherwise we will use the 'RECENTRIFUGE' program to identify contaminants."
            "It is preferred that you align your sequences using the kraken2 'standard' infex file"
    )

    parser_finalisation.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help="The number of cpu threads to use"
    )

    parser_finalisation.set_defaults(func=finalisation)
    
    return parser.parse_args()


def main():
    args = parse_args()
    if len(args.__dict__) <= 1:
        # No arguments or subcommands were given.
        exit("No subcommands provided. Type --help to see what is available")

    args.func(args)

if __name__ == "__main__":
    main()
