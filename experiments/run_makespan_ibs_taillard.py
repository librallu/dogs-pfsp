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

BRANCHING = ["forward", "bi_min"]
GUIDES = ["bound", "idle", "alpha", "wfrontalpha", "gap"]

def create_experiment(n,m,i,branching,guide):
    return lambda: execute_program(output_dir, "../target/release/dogs-pfsp {} {} {} {} {} {}".format(
        "../insts/Taillard/tai{}_{}_{}.txt".format(n,m,i), # inst
        n*m*0.045, # T
        "makespan",
        guide,
        branching,
        output_dir+"/IBS_{}_{}_{}_{}_{}".format(branching,guide,n,m,i) # output files
        ), "log_IBS_{}_{}_{}_{}_{}".format(branching,guide,n,m,i), "completed_list.txt")

if __name__ == "__main__":
    exp_id = get_git_id()
    output_dir = "IBS_makespan_tai__{}".format(exp_id)
    if os.path.exists(output_dir):
        print("DIRECTORY {} ALREADY EXISTS".format(output_dir))
    else:
        print("CREATING RESULT DIRECTORY ({})".format(output_dir))
        os.mkdir(output_dir)
    print("RUNNING EXPERIMENTS")
    tasks = []
    for n,m in TAI_TYPES:
        for g in GUIDES:
            for b in BRANCHING:
                for i in range(0,10):
                    tasks.append(create_experiment(n,m,i,b,g))
    run_parallel(tasks)
    print("=== EXPERIMENT FINISHED: {} ===".format(exp_id))