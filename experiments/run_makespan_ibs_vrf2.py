#!/usr/bin/python3
from experiments import get_git_id, run_parallel, execute_program
import os


N = [700]
M = [60]
I = range(1,4)

BRANCHING = ["bi_min"]
GUIDES = ["gap"]

def create_experiment(n,m,i,branching,guide):
    return lambda: execute_program(output_dir, "../target/release/dogs-pfsp {} {} {} {} {} {}".format(
        "../insts/Large/VFR{}_{}_{}_Gap.txt".format(n,m,i), # inst
        n*m*0.045, # T
        "makespan",
        guide,
        branching,
        output_dir+"/IBS_VFR_{}_{}_{}_{}_{}".format(branching,guide,n,m,i) # output files
        ), "log_IBS_VFR_{}_{}_{}_{}_{}".format(branching,guide,n,m,i), "completed_list.txt")

if __name__ == "__main__":
    exp_id = get_git_id()
    output_dir = "IBS_makespan_VFR2__{}".format(exp_id)
    if os.path.exists(output_dir):
        print("DIRECTORY {} ALREADY EXISTS".format(output_dir))
    else:
        print("CREATING RESULT DIRECTORY ({})".format(output_dir))
        os.mkdir(output_dir)
    print("RUNNING EXPERIMENTS")
    tasks = []
    for n in N:
        for m in M:
            for g in GUIDES:
                for b in BRANCHING:
                    for i in I:
                        tasks.append(create_experiment(n,m,i,b,g))
    run_parallel(tasks)
    print("=== EXPERIMENT FINISHED: {} ===".format(exp_id))