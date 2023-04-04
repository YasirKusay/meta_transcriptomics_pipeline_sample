from libc.stdio cimport FILE, fopen, fgetc, fclose, feof

cdef extern from "stdio.h":
    ssize_t getline(char **, size_t *, FILE *)

def countNumSeqs(str file, int isFasta):

    cdef FILE *fp
    cdef char c
    cdef int count

    count = 0

    fp = fopen(file.encode(), "r")
    while (feof(fp) is False):
        c = fgetc(fp)
        if (isFasta == 1 and c == '>') or ((isFasta != 1) and (c == '@')):
            count = count + 1;

    fclose(fp)
    return count

def countNumLines(str file):

    cdef FILE *fp
    cdef char c
    cdef int count
    cdef int charsOnCurrentLine

    count = 0

    fp = fopen(file.encode(), "r")
    while (feof(fp) is False):
        c = fgetc(fp)
        if (c == "\n"):
            charsOnCurrentLine = 0;
            count = count + 1;
        else:
            charsOnCurrentLine += 1;

    fclose(fp)

    if charsOnCurrentLine > 0:
        count = count + 1;

    return count