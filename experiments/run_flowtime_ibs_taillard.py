#!/usr/bin/python3
from experiments import get_git_id
import os

BRANCHINGS = [
    ["f_flowtime", "alpha"]
]

output_dir_prefix = "IBS_f_flowtime_tai"

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


def create_experiment(n,m,i,branching,guide,output_dir):
    output_prefix = "{}/TA_{}_{}_{}_{}_{}".format(output_dir, n,m,i,branching,guide)
    exec_string = "tsp ../target/release/dogs-pfsp -i ../insts/Taillard/tai{}_{}_{}.txt -s {} -p {} -t {} {} -g {}".format(
        n,m,i, output_prefix+".sol", output_prefix+".perf.json", n*m*0.045, branching, guide
    )
    os.system(exec_string)

if __name__ == "__main__":
    exp_id = get_git_id()
    output_dir = "{}__{}".format(output_dir_prefix, exp_id)
    if os.path.exists(output_dir):
        print("DIRECTORY {} ALREADY EXISTS".format(output_dir))
    else:
        print("CREATING RESULT DIRECTORY ({})".format(output_dir))
        os.mkdir(output_dir)
    for n,m in TAI_TYPES:
        for i in range(0,10):
            for b,g in BRANCHINGS:
                create_experiment(n,m,i,b,g,output_dir)