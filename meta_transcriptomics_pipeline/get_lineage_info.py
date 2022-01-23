from ete3 import NCBITaxa

def get_taxids(input_file):
    taxids = {}
    with open(input_file, "r") as f:
        for line in f:
            curr = line.split()
            taxid = curr[0]
            abundance = curr[1]
            taxids[taxid] = abundance

    return taxids

def get_lineage_info(input_file, output_file):
    taxids = get_taxids(input_file)
    ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    wf = open(output_file, "w")
    ncbi = NCBITaxa()

    for taxid in taxids.keys():
        lineage = ncbi.get_lineage(taxid)
        lineage2ranks = ncbi.get_rank(lineage)
        ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
        wf.write(taxid)
        for rank in ranks:
            r = ranks2lineage.get(rank, 'Unknown')
            wf.write("\t" + r)

        wf.write("\n")

    wf.close()
