# DOGS-PFSP

DOGS implementation for the Permutation FlowShop Problem (makespan and flowtime minimization).

This solver is built using the [DOGS framework](https://github.com/librallu/dogs)

# building & running

this solver requires the rust programming language installed (see [these instructions](https://www.rust-lang.org/learn/get-started) to install it)

## compilation 

`cargo build --release`

## solver help

`./target/release/dogs-pfsp --help`

## usage example (makespan minimization)

### run solver

`./target/release/dogs-pfsp --instance insts/Large/VFR800_60_10_Gap.txt --time 30 --solution tmp/VFR800_60_10.sol fb_makespan --guide walpha --branching bi_min`

 should produce the following output:
 ```
======================================
	 printing solutions in: tmp/VFR800_60_10.sol
 time (s)      nb nodes       objective      Iter           
------------------------------------------------------------
  0.000                                      start D=1      
  0.126        320.4 k        50.086                        
  0.128        320.4 k        50.086         start D=2      
  0.381        960.4 k        49.141                        
  0.383        960.4 k        49.107                        
  0.385        960.4 k        49.107         start D=4      
  0.898        2.2 M          47.719                        
  0.900        2.2 M          47.719         start D=8      
  1.920        4.8 M          47.186                        
  1.922        4.8 M          47.168                        
  1.924        4.8 M          47.168         start D=16     
  3.969        9.9 M          46.394                        
  3.972        9.9 M          46.394         start D=32     
  8.034        20.1 M         46.174                        
  8.036        20.1 M         46.174         start D=64     
  16.235       40.6 M         45.795                        
  16.237       40.6 M         45.788                        
  16.239       40.6 M         45.788         start D=128    

                nb pruned            234

             nb generated         75.5 M
              nb expanded        164.7 k
               nb trashed         117.0 
     avg branching factor            458
                  nb eval          20.0 
                 nb guide         75.5 M
        time searched (s)         30.000
      generated nodes / s          2.5 M
        solutions created           0.0 
```

### check solution

The checker takes as parameters:
 - the instance
 - the solution
 - the objective (0: makespan, 1: flowtime)

`./checker.exe insts/Large/VFR800_60_10_Gap.txt tmp/VFR800_60_10.sol 0`

should produce:
```
Solution value: 45788
```


## usage example (flowtime minimization)

`./target/release/dogs-pfsp --instance insts/Taillard/tai500_20_9.txt --time 30 --solution tmp/tai500_20_9.sol f_flowtime --guide alpha`

 should produce the following output:

```
======================================
	 printing solutions in: tmp/tai500_20_9.sol
 ============== Forward Flowtime(g=Alpha) ==============

 time (s)      nb nodes       objective      Iter           
------------------------------------------------------------
  0.000                                      start D=1      
  0.007        125.2 k        7.274.402                     
  0.008        125.2 k        7.274.402      start D=2      
  0.023        375.2 k        7.099.737                     
  0.025        375.2 k        7.099.737      start D=4      
  0.060        874.7 k        7.000.026                     
  0.062        874.7 k        7.000.004                     
  0.064        874.7 k        7.000.004      start D=8      
  0.128        1.9 M          6.851.859                     
  0.130        1.9 M          6.851.821                     
  0.131        1.9 M          6.851.821      start D=16     
  0.255        3.9 M          6.802.887                     
  0.257        3.9 M          6.802.884                     
  0.258        3.9 M          6.802.884      start D=32     
  0.507        7.9 M          6.751.376                     
  0.509        7.9 M          6.751.325                     
  0.510        7.9 M          6.751.325      start D=64     
  1.023        15.8 M         6.697.901                     
  1.025        15.8 M         6.697.901      start D=128    
  2.034        31.8 M         6.687.834                     
  2.035        31.8 M         6.687.827                     
  2.037        31.8 M         6.687.826                     
  2.038        31.8 M         6.687.826      start D=256    
  4.072        63.8 M         6.654.700                     
  4.073        63.8 M         6.654.685                     
  4.075        63.8 M         6.654.685      start D=512    
  8.182        127.6 M        6.642.686                     
  8.184        127.6 M        6.642.669                     
  8.186        127.6 M        6.642.669      start D=1024   
  16.595       255.1 M        6.635.267                     
  16.596       255.1 M        6.635.256                     
  16.598       255.1 M        6.635.249                     
  16.600       255.1 M        6.635.240                     
  16.602       255.1 M        6.635.240      start D=2048   
                nb pruned           2025


             nb generated        466.2 M
              nb expanded          1.6 M
               nb trashed          2.0 k
     avg branching factor            287
                  nb eval          44.0 
                 nb guide        466.2 M
        time searched (s)         30.000
      generated nodes / s         15.5 M
        solutions created           0.0 
```

### check solution

The checker takes as parameters:
 - the instance
 - the solution
 - the objective (0: makespan, 1: flowtime)

`./checker.exe insts/Taillard/tai500_20_9.txt tmp/tai500_20_9.sol 1`

should produce:
```
Solution value: 6.63524e+06
```


# Using a randomly generated instance

```./target/release/run_rand_inst --nb_jobs 800 --nb_machines 60 --time 30 --perf tmp/perf.json fb_makespan --guide walpha --branching bi_min```

Should produce the following results:

```
	 printing perfs in: tmp/perf.json

 time (s)      dual      nb nodes       objective      Iter           
----------------------------------------------------------------------
  0.000                                                start D=1      
  0.141        42.958    320.4 k        50.253                        
  0.141        42.958    320.4 k        50.253         start D=2      
  0.419        42.958    960.4 k        47.799                        
  0.420        42.958    960.4 k        47.799         start D=4      
  0.970        42.958    2.2 M          47.581                        
  0.971        42.958    2.2 M          47.569                        
  0.971        42.958    2.2 M          47.569         start D=8      
  2.067        42.958    4.8 M          46.880                        
  2.067        42.958    4.8 M          46.880         start D=16     
  4.249        42.958    9.9 M          46.485                        
  4.249        42.958    9.9 M          46.452                        
  4.249        42.958    9.9 M          46.452         start D=32     
  8.656        43.050    20.1 M         45.999                        
  8.656        43.050    20.1 M         45.999         start D=64     
  18.182       43.050    40.6 M         45.771                        
  18.182       43.050    40.6 M         45.754                        
  18.182       43.050    40.6 M         45.744                        
  18.182       43.050    40.6 M         45.744         start D=128    

               dual bound          43050


                nb pruned            232

             nb generated         68.0 M
              nb expanded        145.0 k
               nb trashed         116.0 
     avg branching factor            468
                  nb eval          22.0 
                 nb guide         68.0 M
        time searched (s)         30.000
      generated nodes / s          2.3 M
        solutions created           0.0 
```