import argparse

def check_fail(command, to_remove):
    if (command.returncode != 0):
        exit_now = True
        print(command.stderr)

    if (exit_now):
        print("fail")
        for file in to_remove:
            os.remove(file)
        
    return False    

def run_pipeline(args: argparse.Namespace):

    ##################### FASTP ########################

    path_command = subprocess.run('which fastp', 
                                    check=True,
                                    shell=True, 
                                    capture_output=True
                                )

    fastp_path = path_command.stdout.decode("utf-8").partition('\n')[0]

    if (path_command.returncode != 0):
        raise Exception("fastp is not installed or initialised")

    fastp_command = fastp_path + 
                    " --in1 " + args.inp1 + 
                    " --in2 " + args.inp2 + 
                    "  --dedup  " + # remove duplicates +
                    " --qualified_quality_phred  " + args.qualified_quality_phred + 
                    " --unqualified_percent_limit " + args.unqualified_percent_limit +
                    " --average_qual " + args.average_qual +
                    " --length_required " + args.length_required +
                    # need to consider adapters, should we give the user a chance to add adatpers?
                    " --thread " + args.threads + 
           
    new_command = subprocess.run(fastp_command, shell=True, capture_output=True)
    if check_fail(new_command, [file1, file2]) is False: return None
    


    #################################### SORTMERNA ############################
    
    path_command = subprocess.run('which sortmerna', 
                                    check=True,
                                    shell=True, 
                                    capture_output=True
                                )

    sortmerna_path = path_command.stdout.decode("utf-8").partition('\n')[0]

    if (path_command.returncode != 0):
        raise Exception("sortmerna is not installed or initialised")

    sortmerna_command = sortmerna_path + 
                    " --ref  " + 
                    " --aligned " + # stores rRNA reads, needs to be deleted afterwards
                    " --other " + # stores non rRNA reads/ creates 2 files, one for forward and reverse strand
                    " --fastx " + # output above 2 in fastx
                    " --reads " + + " --reads " # pe reads, use --reads twice
                    " -a " + args.threads
                    # " -e " evalue threshol, default = 1
                    # below 2 are not necessary to edit, keep them low as they are concerned with rRNA alignment, which we are trying to get rid of anyways
                    " --num-alignments 1 " + # 1 = all alignments reaching E value threshold are reported, 0 = the first alignment passing E-value threshold
                    " --best 1 " + # 1 = all high candidate reference sequences will be searched for alignments
                    
    new_command = subprocess.run(sortmerna_command, shell=True, capture_output=True)
    if check_fail(new_command, [file1, file2]) is False: return None