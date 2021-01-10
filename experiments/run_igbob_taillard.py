#!/usr/bin/python3
from experiments import get_git_id, run_parallel, execute_program
import os


TAI_TYPES = [
    (20, 5),
    (20, 10),
    (20, 20),
    (50, 5),
    (50, 10),
    (50, 20),
    (100, 5),
    (100, 10),
    (100, 20),
    (200, 10),
    (200, 20),
    (500, 20)
]

TIMES = [30, 60, 90]


def create_experiment(n,m,t,i):
    return lambda: execute_program(output_dir, "../litcodes/igbob/igbob.exe {} {} {}".format(
        "../insts/Taillard/tai{}_{}_{}.txt".format(n,m,i), # inst
        t, # T
        output_dir+"/perfprofile_igbob_{}_{}_{}_{}".format(t,n,m,i) # output.json
        ), "log_igbob_{}_{}_{}_{}".format(t,n,m,i), "completed_list.txt")

if __name__ == "__main__":
    exp_id = get_git_id()
    output_dir = "igbob_tai__{}".format(exp_id)
    if os.path.exists(output_dir):
        print("DIRECTORY {} ALREADY EXISTS".format(output_dir))
    else:
        print("CREATING RESULT DIRECTORY ({})".format(output_dir))
        os.mkdir(output_dir)
    print("RUNNING EXPERIMENTS")
    tasks = []
    for n,m in TAI_TYPES:
        for t in TIMES:
            for i in range(0,10):
                tasks.append(create_experiment(n,m,t,i))
    run_parallel(tasks)
    print("=== EXPERIMENT FINISHED: {} ===".format(exp_id))