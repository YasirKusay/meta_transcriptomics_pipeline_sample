import logging
import subprocess
import os
import tempfile

log = logging.getLogger(__name__)

def run_shell_command(command):
    processDetails = subprocess.run(command, shell=True)

    if (processDetails.returncode != 0):
        log.error(processDetails.stderr)
        exit(processDetails.returncode)

def check_command_exists(program_name):
    path_command = subprocess.run('which ' + program_name,
                                        check=True,
                                        shell=True,
                                        capture_output=True
                                )

    if (path_command.returncode != 0):
        log.error(program_name + " is not installed or initialised properly")
        exit(1)

def generate_temp_file(extension, working_dir):
    file = tempfile.NamedTemporaryFile(
        suffix=extension,
        dir=working_dir,
        mode="w",
        delete=False
    )

    file.file.close()

    return file.name

def delete_temp_files(files):
    for file in files:
        os.remove(file)

def check_megahit_success(log): 
    megahit_exit_message = subprocess.check_output('tail -n 1 ' + log, shell=True).decode('utf-8').strip('\n').replace("[","").replace("]","")

    if "Exit code" in megahit_exit_message:
        megahit_exit_status = megahit_exit_message.split(" ")[2]
        if megahit_exit_status.isnumeric() and int(megahit_exit_status) != 0:
            return False
    
    return True