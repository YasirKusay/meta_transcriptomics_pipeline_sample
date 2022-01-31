from libc.stdio cimport FILE, fopen, getline, printf, fclose, fprintf
from libc.string cimport strcmp, strlen
from libc.stdlib cimport malloc, free

cdef extern from "stdio.h":
    ssize_t getline(char **, size_t *, FILE *)

def obtain_relevant_taxids(str accession_file, str mapping_file, str write_file):

    cdef int size, curr_iter, len_accession, len1, i
    cdef char *line
    cdef char *curr_accession 
    cdef int cancel
    cdef char *curr = <char *>malloc(100 * sizeof (char))
    cdef FILE *mf 
    cdef FILE *wf
    cdef FILE *cf
    cdef size_t big_size, small_size
    cdef ssize_t big_ssize, small_ssize

    for n in range(100):
        curr[n] = '\0'

    size = 0
    line = NULL

    mf = fopen(accession_file.encode(), "r")
    while (True):
        small_ssize = getline(&curr_accession, &small_size, mf)
        if (small_ssize < 0): 
            break
        size += 1

    fclose(mf)
    
    mf = fopen(accession_file.encode(), "r")
    small_ssize = getline(&curr_accession, &small_size, mf)
    len_accession = strlen(curr_accession)
    if (curr_accession[len_accession - 1] == "\n"):
        curr_accession[len_accession-1] = '\0';
    wf = fopen(write_file.encode(), "w")

    curr_iter = 0

    cf = fopen(mapping_file.encode(), "r")
    big_ssize = getline(&line, &big_size, cf)
    while True:
        big_ssize = getline(&line, &big_size, cf)
        if big_ssize < 0:
            break
        len1 = strlen(line)
        for i in range(0, len1 - 1):
            if line[i] == '\0' or line[i] == '\t' or line[i] == '\n':
                curr[i] = '\0'
                break
            curr[i] = line[i];
        printf("CURR_LINE: %s\n", curr)
        if (strcmp(curr, curr_accession) == 0):
            printf("GENERIC\n")
            fprintf(wf, line)
            curr_iter += 1
            if (curr_iter == size):
                break
            small_ssize = getline(&curr_accession, &small_size, mf)
            len_accession = strlen(curr_accession)
            if (curr_accession[len_accession - 1] == "\n"):
                curr_accession[len_accession-1] = '\0';
           
        elif (strcmp(curr, curr_accession) > 0):
            cancel = 0
            while (strcmp(curr, curr_accession) > 0):
                printf("SKIPPING %s, %s\n", curr_accession, curr)
                curr_iter += 1
                if (curr_iter == size):
                    cancel = 1
                    break
                small_ssize = getline(&curr_accession, &small_size, mf)
                len_accession = strlen(curr_accession)
                if (curr_accession[len_accession - 1] == "\n"):
                    curr_accession[len_accession-1] = '\0';
                printf("NOW: %s\n", curr_accession)
                if (strcmp(curr, curr_accession) == 0):
                    break

            if cancel == 1:
                break

            printf("STRCMP: %s %s\n", curr_accession, curr_accession)
            printf("%d\n", strcmp(curr, curr_accession))

            if (strcmp(curr, curr_accession) == 0):
                printf("HELLO\n")
                fprintf(wf, line)
                printf("EXTRA\n")
                printf("%s\n", curr)
                printf("%s\n", curr_accession)
                printf("HEHE XD\n")
                curr_iter += 1
                if (curr_iter == size):
                    break
                small_ssize = getline(&curr_accession, &small_size, mf)
                len_accession = strlen(curr_accession)
                if (curr_accession[len_accession - 1] == "\n"):
                    curr_accession[len_accession - 1] = '\0';

    fclose(mf)
    fclose(wf)
    fclose(cf)

