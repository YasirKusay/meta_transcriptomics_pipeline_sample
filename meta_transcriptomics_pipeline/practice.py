mf = open("samp4.txt", "r")
curr_accession = mf.readline().rstrip()
wf = open("write_file", "w")

with open("samp3.txt", "r") as cf:
    for line in cf:
        curr = line.split()
        print(curr[0])
        print("ASD")
        print(curr_accession.rstrip())
        if (curr[0] == curr_accession):
            wf.write(line)
            curr_accession = mf.readline().rstrip()

mf.close()
wf.close()


#diamond - > need to talk about the speed of diamond (need to discuss pro/against of different stuff -> trade offs in options, mem sttings)