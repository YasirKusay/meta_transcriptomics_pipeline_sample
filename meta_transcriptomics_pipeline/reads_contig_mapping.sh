# from the megahit program, run this script so that we can get reads mapped to contigs
bbmap/bbwrap.sh ref=final.contigs.fa in=Blank1.fastq in2=Blank2.fastq out=aln.sam.gz kfilter=22 subfilter=15 maxindel=80

# firstly lets look into the file mapping reads to contigs (sam file)
# gets only reads that mapped against contigs
# format: reads contigs
awk '$3 != "*" && $1 !~ /^@/ {print $1 "\t" $3}' bbwrap_sam.sam > reads_contigs_mapping.txt

# lets remove contigs that failed to map entirely
# TO DO, compare final_contigs.fa with the final mapping file





