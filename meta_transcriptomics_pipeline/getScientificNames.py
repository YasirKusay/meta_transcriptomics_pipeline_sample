def getScientificNames(taxids, folder_location):
    ranked_lineage = folder_location + "/names.dmp"
    sci_names = {}

    taxids = list(dict.fromkeys(taxids))
    taxids.sort()
    inc = 0
    found = 0

    # it appears that there are values being skipped, or eof gets reached earlier
    big_inc = 1
    with open(ranked_lineage, "r") as f:
        for line in f:
            curr = line.split("\t|\t")
            if int(curr[0]) == int(taxids[inc]):
                if curr[3] == "scientific name\t|\n":
                    sci_names[curr[0]] = curr[1]
                    inc += 1
                    found += 1
                    if inc == len(taxids):
                        break
            elif int(curr[0]) > int(taxids[inc]): # has been a problem before, e.g. #115782 is not in names.dmp for some reason, simply skips that taxid
                sci_names[curr[0]] = taxids[inc] 
                # if taxid is not in names.dmp, simply have the "scientific name" be its taxid
                inc += 1

                # good idea to check anyways, in case we may accidentally skip it
                
                if int(curr[0]) == int(taxids[inc]):
                    if curr[3] == "scientific name\t|\n":
                        sci_names[curr[0]] = curr[1]
                        inc += 1
                        if inc == len(taxids):
                            break


    if (inc < len(taxids)):
        while (inc < len(taxids)):
            sci_names[taxids[inc]] = taxids[inc]
            inc += 1    


    return sci_names
