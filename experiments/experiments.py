#!/usr/bin/python3
"""
Runs experiments.
 1. creates numerical experiment output directory (date and git hash)
 2. runs all testings (possibly in parallel)
"""

import subprocess
from datetime import date
import multiprocessing
from time import sleep
from concurrent.futures import ThreadPoolExecutor
from subprocess import Popen, PIPE
import os
import shlex
import sys


def get_git_id():
    """
    returns git last commit version (used to track which code was run)
    """
    res = subprocess.check_output(["git", "rev-parse", "HEAD"]).strip()
    return str(res.decode("utf-8"))


def get_time_id():
    """
    returns the current date (used to track when a code was run)
    """
    return str(date.today().strftime("%Y_%m_%d"))


def get_experiment_id():
    """
    returns the experiment id (date+git_id)
    """
    return get_time_id()+"__"+get_git_id()


def run_parallel(jobs):
    """
    run parallel jobs on all available CPUs
    :jobs: list of functions to be called
    """
    nbcpu = multiprocessing.cpu_count()/4 # gets number of physical cpus and divide by 2 to get more accurate results
    print("using {} CPUs".format(nbcpu))
    with ThreadPoolExecutor(max_workers=nbcpu) as executor:
        job_list = [ executor.submit(j) for j in jobs ]
        while True:
            all_completed = True
            for j in job_list:
                if j.running():
                    all_completed = False
            if all_completed:
                break
        sleep(1)  # spend most of the time sleeping to save processor usage
    print("all jobs finished")
        


def execute_program(output_dir, command, logfilename, completedfilename):
    """
    :output_dir: output directory to write files
    :command: command to run (string)
    :logfilename: prints the stdout in the logfilename file
    :completedfilename: dictionnary containing all experiments previously done
    """
    # check if experiment was already done
    completedlist = []
    # if completedfilename does not exists, create it
    if not os.path.exists(output_dir+"/"+completedfilename):
        with open(output_dir+"/"+completedfilename, "w") as f:
            f.write("")
    with open(output_dir+"/"+completedfilename, "r") as f:
        completedlist = f.read().splitlines()
    if command in completedlist:
        sys.stdout.write("EXPERIMENT {} was already done\n".format(command))
    else:
        iscompleted = False
        sys.stdout.write("STARTING EXPERIMENT {}\n".format(command))
        process = Popen(shlex.split(command), stdout=PIPE)
        (out,err) = process.communicate()
        code = process.wait()
        if code == 0:
            iscompleted = True
        # write log into logfilename
        with open(output_dir+"/"+logfilename, "w") as f:
            f.write("##### LOG STDOUT:\n{}\n##### LOG STDERR:\n{}\n".format(out,err))
        # add experiment to the completed dictionnary file 
        if iscompleted:
            with open(output_dir+"/"+completedfilename, "a") as f:
                f.write("{}\n".format(command))
        sys.stdout.write("FINISHED EXPERIMENT {}\n".format(command))



def test_function(i):
    print("START {}".format(i))
    sleep(i)
    print("END {}".format(i))



def main():
    print(get_experiment_id())
    run_parallel([
        lambda: test_function(1),
        lambda: test_function(2),
        lambda: test_function(3),
    ])


if __name__ == "__main__":
    main()