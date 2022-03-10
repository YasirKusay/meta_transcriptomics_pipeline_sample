import os

def filter_files(infile, out1, out2):
    if (os.path.splitext(infile)[1] == ".fq"):
        filter_fq(infile, out1, out2)

def filter_fq(infile, out1, out2):
    line_inc = 1

    f1 = open(out1, "w")
    f2 = open(out2, "w")

    with open(infile) as f:
        curr_name = ""
        curr_seq = ""
        curr_extra = ""
        for line in f:
            if line_inc == 1:
                curr_name = line
            elif line_inc == 2:
                curr_seq = line
            elif line_inc == 4:
                curr_extra = line
                line_inc = 1

                if (len(line.strip()) < 100):
                    f2.write(curr_name)
                    f2.write(curr_seq)
                    f2.write("+\n")
                    f2.write(curr_extra)
                else:
                    f1.write(curr_name)
                    f1.write(curr_seq)
                    f1.write("+\n")
                    f1.write(curr_extra)

                continue

            line_inc += 1
