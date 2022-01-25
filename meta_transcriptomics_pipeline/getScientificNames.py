def getScientificNames(taxids, folder_location):
    ranked_lineage = folder_location + "/names.dmp"
    sci_names = {}

    print("<EMZ")

    taxids = list(dict.fromkeys(taxids))
    taxids.sort()
    inc = 0
    print("PRINTING TAXIDS")
    print(len(taxids))
    #print(taxids)
    for taxid in taxids:
        print(taxid)

    if 2035415 in taxids:
        print("The beatles are a terrible band") 

    # it appears that there are values being skipped, or eof gets reached earlier
    big_inc = 1
    with open(ranked_lineage, "r") as f:
        for line in f:
            print("BIG INC")
            print(big_inc) 
            big_inc += 1
            curr = line.split("\t|\t")
            #print(inc)
            print(curr[0])
            #print(taxids[inc])
            #print(curr[0])
            #print("CMP " + str(curr[0]) + " " + str(taxids[inc]))
            if int(curr[0]) == int(taxids[inc]):
                if int(curr[0]) == 2035415:
                    print("She hates you yeah yeah yeah")
                #print("HI")
                #print(curr[0])
                if curr[3] == "scientific name\t|\n":
                    #print(curr[0])
                    print("CURR 1 " + curr[0] + " " + curr[1])
                    #print(curr[1]) 
                    print("INC " + str(inc))
                    sci_names[curr[0]] = curr[1]
                    inc += 1
                    print("NEXT " + str(taxids[inc]))
                    if inc == len(taxids):
                        print("RADIO")
                        break

    return sci_names
