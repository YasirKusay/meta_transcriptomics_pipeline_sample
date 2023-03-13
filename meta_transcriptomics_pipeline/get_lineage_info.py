from ete3 import NCBITaxa
from meta_transcriptomics_pipeline.getScientificNames import getScientificNames

# returns a dictionary
def get_lineage_info(taxids, taxdump_location):
    ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    ncbi = NCBITaxa()
    all_species_lineages = []
    all_taxids = []

    for curr_taxid in taxids:
        # TODO, NEED A WAY TO DEAL WITH "UNKNOWN"
        curr_lineage_full = []
        if curr_taxid == "Unknown":
            continue
        try:
            curr_lineage_full = ncbi.get_lineage(curr_taxid) # gets lineage of current taxid, returns only the its taxids
        except:
            print(str(curr_taxid) + " NOT FOUND")
            continue
        lineage2ranks = ncbi.get_rank(curr_lineage_full) # gets ranks of the lineages in curr_lineage_full as a dict where rank is the key
        ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items()) # same as above, flips keys however
        curr_lineage = []

        if curr_taxid != "Unknown":
            all_taxids.append(int(curr_taxid))
        for rank in ranks:
            r = ranks2lineage.get(rank, 'Unknown') # returns value of rank if it exists, 'Unknow' otherwise
            if r != "Unknown":
                all_taxids.append(int(r))
            curr_lineage.append(r)
        
        # previously, curr_lineage used to start with the abundance (read/tmp),
        # now it starts with the taxid of the highest level thing of the species
        all_species_lineages.append(curr_lineage) 
    
    sci_names = getScientificNames(all_taxids, taxdump_location) 

    all_lineages_resolved = {}

    for curr_lineage in all_species_lineages:
        curr_species_taxid = curr_lineage[-1]
        all_lineages_resolved[curr_species_taxid] = []

        for item in curr_lineage:
            all_lineages_resolved[curr_species_taxid].append(str(sci_names[str(item)]))

    return all_lineages_resolved
