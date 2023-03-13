def getScientificNames(taxids, taxdump_location):
    ranked_lineage = taxdump_location + "/names.dmp"
    sci_names = {}

    taxids = list(dict.fromkeys(taxids))
    taxids.sort()
    inc = 0

    # it appears that there are values being skipped, or eof gets reached earlier
    # for the file, each line is split by a "|"" where the first column is the taxid
    # the second column is a scientific name
    # the third column I am not sure about (sometimes it is empty)
    # the fourth column contains the information we are after
    # annoyingly, names.dmp has a lot of white space between the columns

    with open(ranked_lineage, "r") as f:
        for line in f:
            curr = line.split("|")
            if int(curr[0].strip()) == int(taxids[inc]): # names.dmp has a lot of whitespace between the | divisions
                if curr[3].strip() == "scientific name":
                    sci_names[str(curr[0].strip())] = curr[1].strip()
                    inc += 1
                    if inc == len(taxids):
                        break
                    if (taxids[inc] == "Unknown"):
                        inc += 1
                    if inc == len(taxids):
                        break
            elif int(curr[0].strip()) > int(taxids[inc]):
                # if the taxid in curr[0] is lexigraphically greater than the current taxid, 
                # it means that it is missing from names.dmp
                # has been a problem before, e.g. #115782 is not in names.dmp for some reason, simply skips that taxid
                # sci_names[str(taxids[inc])] = str(taxids[inc])
                # if taxid is not in names.dmp, simply have the "scientific name" be its taxid
                #inc += 1

                cancel = False
                # need to check if the next few taxids are missing.
                while int(curr[0].strip()) > int(taxids[inc]):
                    sci_names[str(taxids[inc])] = str(taxids[inc])
                    inc += 1
                    if inc == len(taxids):
                        cancel = True
                        break
                    if (taxids[inc] == "Unknown"):
                        inc += 1
                    if inc == len(taxids):
                        cancel = True
                        break

                if (cancel == True):
                    break

                # good idea to check anyways, in case we may accidentally skip it
                
                if int(curr[0].strip()) == int(taxids[inc]):
                    if curr[3].strip() == "scientific name":
                        sci_names[str(curr[0].strip())] = curr[1].strip()
                        inc += 1
                        if inc == len(taxids):
                            break
                        if (taxids[inc] == "Unknown"):
                            inc += 1
                        if inc == len(taxids):
                            break


    if (inc < len(taxids)):
        while (inc < len(taxids)):
            sci_names[str(taxids[inc])] = str(taxids[inc])
            inc += 1    

    sci_names["Unknown"] = "Unknown"

    return sci_names
