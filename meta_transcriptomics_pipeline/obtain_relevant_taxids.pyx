from libc.stdio cimport FILE, fopen, getline, fclose, fprintf
from libc.string cimport strcmp, strlen
# from libcpp cimport bool

cdef extern from "stdio.h":
    #FILE * fopen ( const char * filename, const char * mode )
    #FILE *fopen(const char *, const char *)
    #int fclose ( FILE * stream )
    #int fclose(FILE *)
    #ssize_t getline(char **lineptr, size_t *n, FILE *stream);
    ssize_t getline(char **, size_t *, FILE *)

def obtain_relevant_taxids(str accession_file, str mapping_file, str write_file):

    cdef int size, curr_iter, len_accession, len1, i
    cdef char *line
    cdef char *curr_accession 
    cdef int cancel
    cdef char *curr
    cdef FILE *mf 
    cdef FILE *wf
    cdef FILE *cf
    cdef size_t big_size, small_size
    cdef ssize_t big_ssize, small_ssize

    size = 0

    with open(accession_file, "r") as f:
        for line in f:
            size += 1

    mf = fopen(accession_file, "r")
    # curr_accession = mf.readline().strip()
    small_ssize = getline(&curr_accession, &small_size, mf)
    len_accession = strlen(curr_accession)
    if (curr_accession[len_accession - 1] == "\n"):
        curr_accession[len_accession-1] = '\0';
    #curr_accession = curr_accession[:-1]
    wf = fopen(write_file, "w")

    curr_iter = 0

    #while True:
    #    if (curr_accession is None or curr_accession == "" or curr_accession == "*"):
    #        #curr_accession = mf.readline().rstrip()
    #        small_ssize = getline(&curr_accession, &small_size, mf)
    #        curr_accession = curr_accession[:-1]
    #        curr_iter += 1
    #        if (curr_iter == size):
    #            break
    #    else:
    #        break

    cf = fopen(mapping_file, "r")
    big_ssize = getline(&line, &big_size, cf)
    while True:
        big_ssize = getline(&line, &big_size, cf)
        if big_ssize == -1: 
            break
        len1 = strlen(line)
        for i in range(0, len1 - 1):
            if line[i] == '\0' or line[i] == '\t' or line[i] == '\n':
                line[i] = '\0'
                break
            curr[i] = line[i];
        if (strcmp(curr, curr_accession) == 0):
            #wf.write(line)
            fprintf(wf, line)
            curr_iter += 1
            if (curr_iter == size):
                break
            small_ssize = getline(&curr_accession, &small_size, mf)
            len_accession = strlen(curr_accession)
            if (curr_accession[len_accession - 1] == "\n"):
                curr_accession[len_accession-1] = '\0';
            #curr_accession = curr_accession[:-1]
            
            #if (curr_accession[len]

        #elif min(curr[0], curr_accession) == curr_accession:
        elif (strcmp(curr, curr_accession) > 0):
            cancel = 0
            while (strcmp(curr, curr_accession) > 0):
                curr_iter += 1
                if (curr_iter == size):
                    cancel = 1
                    break
                small_ssize = getline(&curr_accession, &small_size, mf)
                len_accession = strlen(curr_accession)
                if (curr_accession[len_accession - 1] == "\n"):
                    curr_accession[len_accession-1] = '\0';
                #curr_accession = curr_accession[:-1]
                if (strcmp(curr, curr_accession) == 0):
                    break

            if cancel == 1:
                break

            if (strcmp(curr, curr_accession) == 0):
                #wf.write(line)
                fprintf(wf, line)
                curr_iter += 1
                if (curr_iter == size):
                    break
                small_ssize = getline(&curr_accession, &small_size, mf)
                len_accession = strlen(curr_accession)
                if (curr_accession[len_accession - 1] == "\n"):
                    curr_accession[len_accession - 1] = '\0';
                #curr_accession = curr_accession[:-1]

    fclose(mf)
    fclose(wf)
    fclose(cf)
