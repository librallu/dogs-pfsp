# DOGS-PFSP

DOGS implementation for the Permutation FlowShop Problem (makespan and flowtime minimization).

This solver is built using the [DOGS framework](https://github.com/librallu/dogs)

# building & running

## compilation 

`cargo build --release`

## program help

`./target/release/dogs-pfsp --help`

## running solver

`./target/release/dogs-pfsp insts/Large/VFR800_60_10_Gap.txt 30 makespan wfrontalpha bi_min tmp/test `

 should produce the following output:
 ```
  time (s)      nb nodes       objective      IBS            
------------------------------------------------------------
  0.001                                      start D=1      
  0.638        320.4 k        50.086                        
  0.640        320.4 k        50.086         start D=2      
  1.936        960.4 k        49.141                        
  1.938        960.4 k        49.107                        
  1.940        960.4 k        49.107         start D=4      
  4.464        2.2 M          47.719                        
  4.466        2.2 M          47.719         start D=8      
  9.689        4.8 M          47.186                        
  9.692        4.8 M          47.168                        
  9.724        4.8 M          47.168         start D=16     
  19.901       9.9 M          46.394                        
  19.903       9.9 M          46.394         start D=32     

                nb pruned             48

             nb generated         15.2 M
              nb expanded         32.5 k
               nb trashed          24.0 
     avg branching factor            466
                  nb eval           7.0 
                 nb guide         15.2 M
        time searched (s)         30.000
      generated nodes / s        505.2 k
        solutions created           0.0 

```

## checking solution

we provide a helper program to check if the solutions are correct.

`./checker.exe insts/Large/VFR800_60_10_Gap.txt tmp/test.sol.txt 0`

should produce the following output:
```
Solution value: 46394
```
