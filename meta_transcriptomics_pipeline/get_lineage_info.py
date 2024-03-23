import logging
from ete3 import NCBITaxa
from meta_transcriptomics_pipeline.getScientificNames import getScientificNames

log = logging.getLogger(__name__)

# returns a dictionary
def get_lineage_info(taxids, taxdump_location, addSubspeciesLineage=False):
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
            log.warning("The lineage for the taxid: " + str(curr_taxid) + " was not found. Ignoring the taxid.")
            continue

        all_taxids += [str(item) for item in curr_lineage_full]

        # a situation that we must consider is that some taxids may have been merged to another taxid
        # if addSubspeciesLineage is False, this will not be an issue because we will simply use our
        # taxid and overwrite the merged_taxid

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

        # lets skip, as we would still like to obtain the lineages
        if 'species' not in ranks2lineage.keys():
            bad_taxids.append(curr_taxid)
            continue

        # to identify if a species is ranked 

        # taxids in ranks2lineage are ints, we should convert to strings

        skipSeq = False

        if curr_taxid != "Unknown":
            all_taxids.append(str(curr_taxid))
        for rank in ranks:
            rank_taxid = ranks2lineage.get(rank, 'Unknown') # returns value of rank if it exists as a key in the dict, 'Unknown' otherwise
            if rank_taxid != "Unknown":
                all_taxids.append(str(rank_taxid))
            if rank_taxid == "Unknown" and rank == "species":
                bad_taxids.append(curr_taxid)
                skipSeq = True
                break

            # taxid of a hit is not a "species taxid"
            # this is very important to filter out as it could mess up
            # the krona chart filtering step
            # ie, if there is something that is higher than the species rank and is selected to 
            # be filtered out, it could remove anything under it, even if they pass the filtering 
            # step, this was the case for one of our samples
            if rank_taxid == curr_taxid and rank != "species":
                bad_taxids.append(curr_taxid)
                skipSeq = True
                break

            # sometimes our taxid's rank is not a species (similar to the scenario above),
            # however, it does not necessarily mean that its rank is a species
            # it can be a subspecies of the species (i.e. lower rank)
            # in this case, prioritise the taxids rank over the actual 'species' rank
            if rank != "species":
                curr_lineage.append(str(rank_taxid))
            else:
                if rank_taxid == curr_taxid:
                    curr_lineage.append(str(curr_taxid))
                else:
                    if addSubspeciesLineage is False:
                        curr_lineage.append(str(curr_taxid)) 
                    else:
                        # we would like add the species lineage and anything after that
                        # go through curr_lineage_full, get indexes of rank_taxid and full_taxid

                        try:
                            for j in range(curr_lineage_full.index(rank_taxid), curr_lineage_full.index(int(curr_taxid))):
                                curr_lineage.append(str(curr_lineage_full[j])) 
                        except:
                            pass
                        curr_lineage.append(str(curr_taxid))

        # previously, curr_lineage used to start with the abundance (read/tmp),
        # now it starts with the taxid of the highest level thing of the species
        if skipSeq is False:
            all_species_lineages.append(curr_lineage) 
    
    # simply add the full lineage
    if addSubspeciesLineage is True:
        for curr_bad_taxid in bad_taxids:
            curr_lineage_full = []
            if curr_bad_taxid == "Unknown":
                continue
            try:
                curr_lineage_full = ncbi.get_lineage(curr_bad_taxid) # gets lineage of current taxid, returns only the its taxids
            except:
                continue
            all_species_lineages.append([str(item) for item in curr_lineage_full])

    # all_taxids will NEVER contain 'Unknown'
    sci_names = getScientificNames(all_taxids, taxdump_location) 

    all_lineages_resolved = {}

    for curr_lineage in all_species_lineages:
        curr_species_taxid = curr_lineage[-1] # curr_species_taxid is a string
        all_lineages_resolved[curr_species_taxid] = []

        for item in curr_lineage:
            all_lineages_resolved[curr_species_taxid].append(sci_names[item].replace(",", " ")) # comma could mess up the Krona chart later, as the scores are comma delimited

    return all_lineages_resolved, bad_taxids