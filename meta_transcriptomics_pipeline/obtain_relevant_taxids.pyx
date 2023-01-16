from libc.stdio cimport FILE, fopen, getline, fclose, fprintf
from libc.string cimport strcmp, strlen, strncpy
from libc.stdlib cimport malloc, free

cdef extern from "stdio.h":
    ssize_t getline(char **, size_t *, FILE *)

def obtain_relevant_taxids(str accession_file, str mapping_file, str write_file):

    cdef int curr_iter, len_accession, len1, i
    cdef char *curr_mapping_file_line = <char *>malloc(100 * sizeof (char))
    cdef char *curr_accession = <char *>malloc(100 * sizeof (char))
    cdef int cancel
    cdef char *curr_mapping_file_accession = <char *>malloc(100 * sizeof (char))
    cdef FILE *af 
    cdef FILE *wf
    cdef FILE *mf
    cdef size_t big_size, small_size
    cdef ssize_t mapping_file_bytes, accession_file_bytes

    for n in range(100):
        curr_mapping_file_accession[n] = '\0'

    size = 0
    curr_mapping_file_line = NULL

    curr_accession = NULL
    small_size = 0

    af = fopen(accession_file.encode(), "r")
    accession_file_bytes = getline(&curr_accession, &small_size, af)
    if (accession_file_bytes < 0):
        return
    len_accession = strlen(curr_accession)
    if (curr_accession[len_accession - 1] == "\n"):
        curr_accession[len_accession - 1] = '\0';
    
    wf = fopen(write_file.encode(), "w")

    curr_iter = 0

    mf = fopen(mapping_file.encode(), "r")
    curr_mapping_file_line = NULL
    big_size = 0
    mapping_file_bytes = getline(&curr_mapping_file_line, &big_size, mf)

    while True:
        mapping_file_bytes = getline(&curr_mapping_file_line, &big_size, mf)
        if mapping_file_bytes < 0:
            break
        len1 = strlen(curr_mapping_file_line)
        for i in range(0, len1 - 1):
            if curr_mapping_file_line[i] == '\0' or curr_mapping_file_line[i] == '\n':
                break
            if curr_mapping_file_line[i] == '\t':
                strncpy(curr_mapping_file_accession, curr_mapping_file_line, i);
                curr_mapping_file_accession[i] = "\0"
                break
        if (strcmp(curr_mapping_file_accession, curr_accession) == 0):
            fprintf(wf, curr_mapping_file_line)
            accession_file_bytes = getline(&curr_accession, &small_size, af)
            if (accession_file_bytes < 0):
                break
            len_accession = strlen(curr_accession)
            if (curr_accession[len_accession - 1] == "\n"):
                curr_accession[len_accession - 1] = '\0';
           
        elif (strcmp(curr_mapping_file_accession, curr_accession) > 0):
            cancel = 0
            while (strcmp(curr_mapping_file_accession, curr_accession) > 0):
                accession_file_bytes = getline(&curr_accession, &small_size, af)
                if (accession_file_bytes < 0):
                    cancel = 1
                    break
                len_accession = strlen(curr_accession)
                if (curr_accession[len_accession - 1] == "\n"):
                    curr_accession[len_accession - 1] = '\0';
                if (strcmp(curr_mapping_file_accession, curr_accession) == 0):
                    break

            if cancel == 1:
                break

            if (strcmp(curr_mapping_file_accession, curr_accession) == 0):
                fprintf(wf, curr_mapping_file_line)
                accession_file_bytes = getline(&curr_accession, &small_size, af)
                if (accession_file_bytes < 0):
                    break
                len_accession = strlen(curr_accession)
                if (curr_accession[len_accession - 1] == "\n"):
                    curr_accession[len_accession - 1] = '\0';

    fclose(af)
    fclose(wf)
    fclose(mf)

