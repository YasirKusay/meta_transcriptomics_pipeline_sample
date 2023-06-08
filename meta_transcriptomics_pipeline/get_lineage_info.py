from ete3 import NCBITaxa
from meta_transcriptomics_pipeline.getScientificNames import getScientificNames

# returns a dictionary
def get_lineage_info(taxids, taxdump_location):
    ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    ncbi = NCBITaxa()
    all_species_lineages = []
    all_taxids = []
    bad_taxids = []

    for curr_taxid in taxids:
        # TO DO, NEED A WAY TO DEAL WITH "UNKNOWN"
        curr_lineage_full = []
        if curr_taxid == "Unknown":
            continue
        try:
            curr_lineage_full = ncbi.get_lineage(curr_taxid) # gets lineage of current taxid, returns only the its taxids
        except:
            print("The lineage for the taxid: " + str(curr_taxid) + " was not found. Ignoring the taxid.")
            continue

        lineage2ranks_unsorted = ncbi.get_rank(curr_lineage_full)  # gets ranks of the lineages in curr_lineage_full as a dict where taxid is the key
        # sorting dictionary values by the original order in curr_lineage full, lineage2ranks_unsorted keys is ranked in ascending order
        lineage2ranks = {i: lineage2ranks_unsorted[i] for i in curr_lineage_full}
        ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items()) # same as above, flips keys however
        curr_lineage = []

        # need to identify if the taxid of a hit is not ranked 'species'
        # this is very important to filter out as it could mess up
        # the krona chart filtering step
        # ie, if there is something that is higher than the species rank and is selected to 
        # be filtered out, it could remove anything under it, even if they pass the filtering 
        # step, this was the case for one of our samples
        # if the below is true, then our taxid is not a species

        if 'species' not in ranks2lineage.keys():
            bad_taxids.append(curr_taxid)
            continue

        # to identify if a species is ranked 

        # taxids in ranks2lineage are ints, we should convert to strings

        if curr_taxid != "Unknown":
            all_taxids.append(str(curr_taxid))
        for rank in ranks:
            rank_taxid = ranks2lineage.get(rank, 'Unknown') # returns value of rank if it exists as a key in the dict, 'Unknown' otherwise
            if rank_taxid != "Unknown":
                all_taxids.append(str(rank_taxid))
            if rank_taxid == "Unknown" and rank == "species":
                bad_taxids.append(curr_taxid)

            # taxid of a hit is not a "species taxid"
            # this is very important to filter out as it could mess up
            # the krona chart filtering step
            # ie, if there is something that is higher than the species rank and is selected to 
            # be filtered out, it could remove anything under it, even if they pass the filtering 
            # step, this was the case for one of our samples
            if rank_taxid == curr_taxid and rank != "species":
                bad_taxids.append(curr_taxid)

            # sometimes our taxid's rank is not a species (similar to the scenario above),
            # however, it does not necessarily mean that its rank is a species
            # it can be a subspecies of the species (i.e. lower rank)
            # in this case, prioritise the taxids rank over the actual 'species' rank
            if rank != "species":
                curr_lineage.append(str(rank_taxid))
            else:
                curr_lineage.append(str(curr_taxid)) 

        # previously, curr_lineage used to start with the abundance (read/tmp),
        # now it starts with the taxid of the highest level thing of the species
        all_species_lineages.append(curr_lineage) 

    # all_taxids will NEVER contain 'Unknown'
    sci_names = getScientificNames(all_taxids, taxdump_location) 

    all_lineages_resolved = {}

    for curr_lineage in all_species_lineages:
        curr_species_taxid = curr_lineage[-1] # curr_species_taxid is a string
        all_lineages_resolved[curr_species_taxid] = []

        for item in curr_lineage:
            all_lineages_resolved[curr_species_taxid].append(sci_names[item].replace(",", " "))

    return all_lineages_resolved, bad_taxids