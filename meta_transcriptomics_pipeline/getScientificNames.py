def getScientificNames(taxids, folder_location):
    ranked_lineage = folder_location + "/names.dmp"
    sci_names = {}

    for item in taxids:
        print(item)

    taxids = list(dict.fromkeys(taxids))
    taxids.sort()
    inc = 0
    found = 0

    # it appears that there are values being skipped, or eof gets reached earlier
    big_inc = 1
    with open(ranked_lineage, "r") as f:
        for line in f:
            curr = line.split("\t|\t")
            #print(line)
            #print(len(taxids))
            #print(taxids[inc-1])
            if int(curr[0]) == int(taxids[inc]):
                if curr[3] == "scientific name\t|\n":
                    sci_names[str(curr[0])] = curr[1]
                    inc += 1
                    found += 1
                    if inc == len(taxids):
                        break
            elif int(curr[0]) > int(taxids[inc]): # has been a problem before, e.g. #115782 is not in names.dmp for some reason, simply skips that taxid
                #sci_names[str(taxids[inc])] = str(taxids[inc])
                # if taxid is not in names.dmp, simply have the "scientific name" be its taxid
                #inc += 1

                cancel = False
                # need to check if the next few taxids are missing.
                while int(curr[0]) > int(taxids[inc]):
                    sci_names[str(taxids[inc])] = str(taxids[inc])
                    inc += 1
                    if inc == len(taxids):
                        cancel = True
                        break

                if (cancel == True):
                    break

                # good idea to check anyways, in case we may accidentally skip it
                
                if int(curr[0]) == int(taxids[inc]):
                    if curr[3] == "scientific name\t|\n":
                        sci_names[str(curr[0])] = curr[1]
                        inc += 1
                        if inc == len(taxids):
                            break


    if (inc < len(taxids)):
        while (inc < len(taxids)):
            sci_names[str(taxids[inc])] = str(taxids[inc])
            inc += 1    

    for item in sci_names.keys():
        print(str(item) + " " + str(sci_names[item]))

    sci_names["Unknown"] = "Unknown"

    return sci_names
