import subprocess
import os
import tempfile

def check_fail(program_name, command):
    if (command.returncode != 0):
        print("########### COMMMAND FAILED ###########")
        print(program_name + " failed")
        print(command.stderr)
        exit(command.returncode)

def check_command_exists(program_name):
    path_command = subprocess.run('which ' + program_name,
                                        check=True,
                                        shell=True,
                                        capture_output=True
                                )

    if (path_command.returncode != 0):
        print(program_name + " is not installed/initialised properly")
        return False

    return True

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
