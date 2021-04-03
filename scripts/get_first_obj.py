#!/usr/bin/python3
import json
from sys import argv

if __name__ == "__main__":
    filename = argv[1]
    print(filename)
    with open(filename) as f:
        perfprofile = json.load(f)
        print("t:{:.3}\tv:{}".format(perfprofile["points"][0]["t"], perfprofile["points"][0]["v"]))