import os

# separates reads inside the infile by its lenth
# longer (> 100bp) reads go into longer_reads
# the rest go into shorter_reads

def separate_reads_by_size(input_reads, longer_reads, shorter_reads):
    if (os.path.splitext(input_reads)[1] == ".fq"):
        separate_reads_fq(input_reads, longer_reads, shorter_reads)

def separate_reads_fq(input_reads, longer_reads, shorter_reads):
    line_inc = 1 # used to determine what part of the current read we are in (the name, sequence, + or quality)

    f1 = open(longer_reads, "w")
    f2 = open(shorter_reads, "w")

    with open(input_reads) as f:
        curr_read_name = ""
        curr_seq = ""
        curr_quality_score = ""
        for line in f:
            if line_inc == 1:
                curr_read_name = line
            elif line_inc == 2:
                curr_seq = line
            elif line_inc == 4:
                curr_quality_score = line
                line_inc = 1

                if (len(line.strip()) < 100):
                    f2.write(curr_read_name)
                    f2.write(curr_seq)
                    f2.write("+\n")
                    f2.write(curr_quality_score)
                else:
                    f1.write(curr_read_name)
                    f1.write(curr_seq)
                    f1.write("+\n")
                    f1.write(curr_quality_score)

                continue

            line_inc += 1
