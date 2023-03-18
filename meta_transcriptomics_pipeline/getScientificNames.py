def getScientificNames(taxids, taxdump_location):
    ranked_lineage = taxdump_location + "/names.dmp"
    sci_names = {}

    taxids = list(dict.fromkeys(taxids))

    # remember that taxids is a list of strings
    # need to convert list to int to sort
    # then convert it back to string
    # since taxids in file are strings, it will make comparisons much clearer
    tmpList = [int(taxid) for taxid in taxids]
    tmpList.sort()
    taxids = [str(taxid) for taxid in tmpList]

    # list increment of taxids list
    inc = 0

    # it appears that there are values being skipped, or eof gets reached earlier
    # for the file, each line is split by a "|"" where the first column is the taxid
    # the second column is a scientific name
    # the third column I am not sure about (sometimes it is empty)
    # the fourth column contains the information we are after, it describes the type of name in column 2
    # annoyingly, names.dmp has a lot of white space between the columns

    with open(ranked_lineage, "r") as f:
        for line in f:
            curr = line.split("|")
            if curr[0].strip() == taxids[inc]: # names.dmp has a lot of whitespace between the | divisions
                if curr[3].strip() == "scientific name":
                    sci_names[taxids[inc]] = curr[1].strip()
                    inc += 1
                    if inc == len(taxids):
                        break
            elif int(curr[0].strip()) > int(taxids[inc]):
                # if the taxid in curr[0] is lexigraphically greater than the current taxid, 
                # it means that it is missing from names.dmp
                # has been a problem before, e.g. #115782 is not in names.dmp for some reason, simply skips that taxid
                # sci_names[str(taxids[inc])] = str(taxids[inc])
                # if taxid is not in names.dmp, simply have the "scientific name" be its taxid
                # inc += 1

                cancel = False
                # need to check if the next few taxids are missing.
                while int(curr[0].strip()) > int(taxids[inc]):
                    sci_names[taxids[inc]] = taxids[inc]
                    inc += 1
                    if inc == len(taxids):
                        cancel = True
                        break

                if (cancel == True):
                    break

                # good idea to check anyways, in case we may accidentally skip it
                
                if curr[0].strip() == taxids[inc]:
                    if curr[3].strip() == "scientific name":
                        sci_names[taxids[inc]] = curr[1].strip()
                        inc += 1
                        if inc == len(taxids):
                            break

    while (inc < len(taxids)):
        sci_names[taxids[inc]] = taxids[inc]
        inc += 1    

    # with our method, we never need to worry about "Unknown being inputted"
    sci_names["Unknown"] = "Unknown"

    return sci_names
