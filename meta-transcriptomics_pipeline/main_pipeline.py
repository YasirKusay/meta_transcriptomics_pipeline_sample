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

def check_command_exists(program_name):
    path_command = subprocess.run('which ' + program_name,
                                        check=True,
                                        shell=True,
                                        capture_output=True
                                )

    if (path_command.returncode != 0):
        path_command = subprocess.run('module add software ' + program_name,
                                        check=True,
                                        shell=True,
                                        capture_output=True
                                )
        if (path_command.returncode != 0):
            raise Exception(program_name + " is not installed or initialised")

        path_command = subprocess.run('which ' + program_name,
                                        check=True,
                                        shell=True,
                                        capture_output=True
                                )

        return path_command

def run_pipeline(args: argparse.Namespace):

    ##################### FASTP ########################

    path_command = check_command_exists('fastp')
    fastp_path = path_command.stdout.decode("utf-8").partition('\n')[0]

    # --dedup removes dups
    fastp_command = fastp_path +\
                    " --in1 " + args.inp1 +\
                    " --in2 " + args.inp2 +\
                    "  --dedup  " +\
                    " --qualified_quality_phred  " + args.qualified_quality_phred +\
                    " --unqualified_percent_limit " + args.unqualified_percent_limit +\
                    " --average_qual " + args.average_qual +\
                    " --length_required " + args.length_required +\
                    " --thread " + args.threads
                    # need to consider adapters, should we give the user a chance to add adatpers?

    new_command = subprocess.run(fastp_command, shell=True, capture_output=True)
    if check_fail(new_command, [file1, file2]) is False: return None



    #################################### SORTMERNA ############################

    path_command = check_command_exists('sortmerna')
    sortmerna_path = path_command.stdout.decode("utf-8").partition('\n')[0]
                       
    #    " --aligned " + # stores rRNA reads, needs to be deleted afterwards
    #                " --other " + # stores non rRNA reads/ creates 2 files, one for forward and reverse strand
    #                " --fastx " + # output above 2 in fastx
    #                " --reads " + + " --reads " # pe reads, use --reads twice
    #                " -a " + args.threads
    #                # " -e " evalue threshol, default = 1
    #                # below 2 are not necessary to edit, keep them low as they are concerned with rRNA alignment, which we are trying to get rid of anyways
    #                " --num-alignments 1 " + # 1 = all alignments reaching E value threshold are reported, 0 = the first alignment passing E-value threshold
    #                " --best 1 " #+ # 1 = all high candidate reference sequences will be searched for alignments

    sortmerna_command = sortmerna_path +\
                    " --ref  " + args.sortmerna_index\
                    " --aligned " +\
                    " --other " +\
                    " --fastx " +\
                    " --reads " + + " --reads "\
                    " -a " + args.threads +\
                    " --num-alignments 1 " +\
                    " --best 1 " # 1 = all high candidate reference sequences will be searched for alignments

    new_command = subprocess.run(sortmerna_command, shell=True, capture_output=True)
    if check_fail(new_command, [file1, file2]) is False: return None
