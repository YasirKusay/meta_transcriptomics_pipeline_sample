import subprocess
import os
import tempfile

def check_fail(command, to_remove):
    if (command.returncode != 0):
        exit_now = True
        print(command.stderr)

    if (exit_now):
        print(command + " failed")
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

    return file

def delete_temp_files(files):
    for file in files:
        os.remove(file)