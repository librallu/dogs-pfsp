#!/usr/bin/python3
from experiments import get_git_id, run_parallel, execute_program
import os


N = [100,200,300,400,500,600,700,800]
M = [20,40,60]
I = range(1,11)

TIMES = [30, 60, 90]


def create_experiment(n,m,t,i):
    return lambda: execute_program(output_dir, "../litcodes/igbob/igbob.exe {} {} {}".format(
        "../insts/Large/VFR{}_{}_{}_Gap.txt".format(n,m,i), # inst
        t, # T
        output_dir+"/perfprofile_igbob_{}_{}_{}_{}".format(t,n,m,i) # output.json
        ), "log_igbob_{}_{}_{}_{}".format(t,n,m,i), "completed_list.txt")


if __name__ == "__main__":
    exp_id = get_git_id()
    output_dir = "igbob_vrf__{}".format(exp_id)
    if os.path.exists(output_dir):
        print("DIRECTORY {} ALREADY EXISTS".format(output_dir))
    else:
        print("CREATING RESULT DIRECTORY ({})".format(output_dir))
        os.mkdir(output_dir)
    print("RUNNING EXPERIMENTS")
    tasks = []
    for n in N:
        for m in M:
            for t in TIMES:
                for i in I:
                    tasks.append(create_experiment(n,m,t,i))
    run_parallel(tasks)
    print("=== EXPERIMENT FINISHED: {} ===".format(exp_id))