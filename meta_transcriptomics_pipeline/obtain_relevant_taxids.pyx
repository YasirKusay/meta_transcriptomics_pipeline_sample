from libc.stdio cimport FILE, fopen, getline, printf, fclose, fprintf
from libc.string cimport strcmp, strlen
from libc.stdlib cimport malloc, free
#from lib.stdlib cimport malloc, free
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
    cdef char *curr = <char *>malloc(100 * sizeof (char))
    cdef FILE *mf 
    cdef FILE *wf
    cdef FILE *cf
    cdef size_t big_size, small_size
    cdef ssize_t big_ssize, small_ssize

    #curr_accession = <char *> malloc((99 + 1) * sizeof(char))
    #curr = malloc(100 * sizeof (char));
    # this below and <char *> above allocates mem to curr[n], otherwise curr[i] - ... below will not work 
    for n in range(100):
        curr[n] = '\0'

    
    size = 0
    line = NULL

    mf = fopen(accession_file.encode(), "r")
    print("MEME")
    while (True):
        #print("JHO")
        small_ssize = getline(&curr_accession, &small_size, mf)
        if (small_ssize < 0): 
            break
        size += 1

    fclose(mf)

    print("SPEED")

    #with open(accession_file, "r") as f:
        #for line in f:
    #        size += 1

    mf = fopen(accession_file.encode(), "r")
    # curr_accession = mf.readline().strip()
    small_ssize = getline(&curr_accession, &small_size, mf)
    len_accession = strlen(curr_accession)
    if (curr_accession[len_accession - 1] == "\n"):
        curr_accession[len_accession-1] = '\0';
    #curr_accession = curr_accession[:-1]
    wf = fopen(write_file.encode(), "w")

    curr_iter = 0

    print("HEHE XD")
    print(mapping_file.encode())

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

    cf = fopen(mapping_file.encode(), "r")
    if cf == NULL:
        print("ASDASDASDASD")
    big_ssize = getline(&line, &big_size, cf)
    #big_ssize = getline(&line, &big_size, cf)
    while True:
        #print("HEHEZ asasdad")
        big_ssize = getline(&line, &big_size, cf)
        #print("CHECKPOINT 0")
        if big_ssize < 0: 
            break
        #print("CHECKPOINT 1")
        len1 = strlen(line)
        #printf("%s\n", line)
        #print("CJECKPOINT 2")
        for i in range(0, len1 - 1):
            #print("I am speedz")
            #if len1 == 0:
            #print("SAD XD")
            if line[i] == '\0' or line[i] == '\t' or line[i] == '\n':
                curr[i] = '\0'
                break
            #print("Checkpoint 2.3")
            #printf("%s\n", line)
            #printf("%c\n", line[i])
            #print("hez")
            curr[i] = line[i];
            #print("Checkpoint 2.4")
        #print("Checkpoint 3")
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

    print("I am done")

    fclose(mf)
    fclose(wf)
    fclose(cf)

    print("Sad XD")
