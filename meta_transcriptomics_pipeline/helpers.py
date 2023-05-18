import subprocess
import os
import tempfile

def run_shell_command(command):
    processDetails = subprocess.run(command, shell=True)

    if (processDetails.returncode != 0):
        print("########### COMMMAND FAILED ###########")
        print(processDetails.stderr)
        exit(processDetails.returncode)

def check_command_exists(program_name):
    path_command = subprocess.run('which ' + program_name,
                                        check=True,
                                        shell=True,
                                        capture_output=True
                                )

    if (path_command.returncode != 0):
        print(program_name + " is not installed or initialised properly")
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
