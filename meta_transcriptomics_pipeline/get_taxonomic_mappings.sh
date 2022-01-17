# we start with the final mapping file for reads/contigs, then proceed to get the tax stuff for both of them
# for each file, we generate 2 files with the following format: read/contig family  genus   species lineage

# lets firstly get the accession info mapped to the reads/contigs
# below has format, contigid accession
awk '$3 != "*" && $1 !~ /^@/ {print $1 "\t" $3}' final_contigs.sam > contigs_accessions_mapping.txt # keep contigs separated instead, so you can check ncbi nt/nr e.g:
# awk '$3 != "*" && $1 !~ /^@/ {print $1 "\t" $3}' contigs_nt.sam > contigs_accessions_mapping.txt
# awk '$3 != "*" && $1 !~ /^@/ {print $1 "\t" $3}' contigs_nr.sam > contigs_accessions_mapping.txt
awk '$3 != "*" && $1 !~ /^@/ {print $1 "\t" $3}' final_reads.sam > reads_accessions_contigs_mapping.txt

# now map reads against contigs
# will have format readid accession

join -1 2 -2 1 -o "1.1, 2.2" <(sort -k2 reads_contigs_mapping.txt) <(sort -k1 contigs_accessions_mapping_nucl.txt) > reads_mapped_sciencename_nucl.txt
join -1 2 -2 1 -o "1.1, 2.2" <(sort -k2 reads_contigs_mapping.txt) <(sort -k1 contigs_accessions_mapping_prot.txt) > reads_mapped_sciencename_prot.txt

awk '{print $2}' reads_mapped_sciencename_nucl.txt | count -c > reads_count.txt

: <<'END_COMMENT'
Alternative, mapping to taxid immediately to reads/contigs

# map taxid to nucl/protein
# can be obtained from: ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/
# has format: query accession taxid
# files need to be sorted
# needs to be done 2 times for each read/contig program, one for nucl, one for prot mapping files
join -1 2 -2 1 -o "1.1, 1.2, 2.1" <(sort -k2 contigs_accessions_mapping.txt) <(sort -k1 ../../../seqstaxids.map) > contigs_accessions_nucl.txt
join -1 2 -2 1 -o "1.1, 1.2, 2.1" <(sort -k2 reads_accessions_contigs_mapping.txt) <(sort -k1 ../../../seqstaxids.map) > reads_accessions_nucl.txt
join -1 2 -2 1 -o "1.1, 1.2, 2.1" <(sort -k2 contigs_accessions_mapping.txt) <(sort -k1 ../../../seqstaxids_prot.map) > contigs_accessions_prot.txt
join -1 2 -2 1 -o "1.1, 1.2, 2.1" <(sort -k2 reads_accessions_contigs_mapping.txt) <(sort -k1 ../../../seqstaxids_prot.map) > reads_accessions_prot.txt

# now we have taxid, we can now map family, name, lineage, etc
# we can map taxid to scientific name using names.dmp
# has format: query accession taxid sci_name
join -1 3 -2 1 -o "1.1, 1.2, 1.3, 2.4" <(sort -k3 contigs_accessions_nucl.txt) <(awk '$1=$1' FS="|" OFS="\t" /data/bio/ncbi_tax/names.dmp | sort -k1) > contigs_sciencename_nucl.txt
join -1 3 -2 1 -o "1.1, 1.2, 1.3, 2.4" <(sort -k3 reads_accessions_nucl.txt) <(awk '$1=$1' FS="|" OFS="\t" /data/bio/ncbi_tax/names.dmp | sort -k1) > reads_sciencename_nucl.txt
join -1 3 -2 1 -o "1.1, 1.2, 1.3, 2.4" <(sort -k3 contigs_accessions_prot.txt) <(awk '$1=$1' FS="|" OFS="\t" /data/bio/ncbi_tax/names.dmp | sort -k1) > contigs_sciencename_prot.txt
join -1 3 -2 1 -o "1.1, 1.2, 1.3, 2.4" <(sort -k3 reads_accessions_prot.txt) <(awk '$1=$1' FS="|" OFS="\t" /data/bio/ncbi_tax/names.dmp | sort -k1) > reads_sciencename_prot.txt

# now we must get the lineage info
# use nodes.dmp to get parent taxid, then map that taxid using names.dmp
# TO DO LATER, WORRY ABOUT IT AFTER U GET TO THE VISUALISATION STAGE :)

# now we need to map the reads against contigs
# will have format readid taxid sciname
join -1 2 -2 1 -o "1.1, 2.3, 2.4" <(sort -k2 reads_contigs_mapping.txt) <(sort -k1 contigs_sciencename_nucl.txt) > reads_mapped_sciencename_nucl.txt
join -1 2 -2 1 -o "1.1, 2.3, 2.4" <(sort -k2 reads_contigs_mapping.txt) <(sort -k1 contigs_sciencename_prot.txt) > reads_mapped_sciencename_prot.txt

# now need to combine the files into one single file to perform TPM calculations

# INSTEAD OF ABOVE, combine 
END_COMMENT
