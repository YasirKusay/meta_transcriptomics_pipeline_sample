import subprocess
import os
import tempfile

def check_fail(program_name, command, to_remove):
    exit_now = False
    if (command.returncode != 0):
        exit_now = True
        print("########### COMMMAND FAILED ###########")

    if (exit_now):
        print(program_name + " failed")
        for f in to_remove:
            os.remove(f)
        return True

    return False

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
