use std::cmp::max;
use ordered_float::OrderedFloat;
use dogs::searchspace::{SearchSpace, SearchTree, TotalChildrenExpansion, GuidedSpace};
use crate::nehhelper::{compute_eqf, EQF, compute_idle_time};

use crate::pfsp::{JobId, Time, Instance};


pub enum TiebreakingScheme {
    NO,  // no tie breaking
    FF   // FF tie breaking (see https://www.sciencedirect.com/science/article/abs/pii/S0305054813003638) 
}


#[derive(Debug, Clone)]
pub struct NEHNode {
    inserted_jobs: Vec<JobId>, // already inserted jobs (partial sequence)
    bound: Time,
    idletime: Option<i64>,  // optionnal idle time of the node (used in tie breaking for instance)
}

pub struct NEHSearch {
    inst: Instance,
    ordered_jobs: Vec<JobId>,  // ordered jobs according to the ordering heuristic
    tie_breaking: TiebreakingScheme,
}


impl SearchSpace<NEHNode, Vec<JobId>> for NEHSearch {
    fn solution(&mut self, node: &NEHNode) -> Vec<JobId> {
        debug_assert!(self.goal(node));
        return node.inserted_jobs.clone();
    }
}


impl GuidedSpace<NEHNode, OrderedFloat<f64>> for NEHSearch {
    fn guide(&mut self, node: &NEHNode) -> OrderedFloat<f64> {
        // let alpha:f64 = (node.inserted_jobs.len() as f64)/(self.inst.nb_jobs() as f64);
        match self.tie_breaking {
            TiebreakingScheme::NO => OrderedFloat(node.bound as f64),
            TiebreakingScheme::FF => match node.idletime {
                None => OrderedFloat(node.bound as f64),
                Some(v) => OrderedFloat((node.bound as f64)-1./(v as f64))
                // Some(v) => OrderedFloat(
                //     (node.bound as f64)*alpha +
                //     ((v as f64)/(self.inst.nb_machines() as f64))*(1.-alpha)
                // )
            }
        }
    }
}


impl TotalChildrenExpansion<NEHNode> for NEHSearch {
    fn children(&mut self, node: &mut NEHNode) -> Vec<NEHNode> {
        debug_assert!(!self.goal(node));  // should not be goal
        let mut res:Vec<NEHNode> = Vec::new();
        let k = node.inserted_jobs.len();  // kth job to be inserted in the partial solution
        let j = self.ordered_jobs[k];  // job to be inserted
        let m = self.inst.nb_machines();
        if node.inserted_jobs.is_empty() {
            // if first job, just insert the first job in the ordered list
            let mut bound = 0;
            for i in 0..m {
                bound += self.inst.p(j, i);
            }
            res.push(NEHNode { 
                inserted_jobs: vec![j],
                bound: bound,
                idletime: None,
            })
        } else {
            let mut eqf:EQF = compute_eqf(&self.inst, &node.inserted_jobs, j);
            for i in 0..m {
                eqf.q[i as usize].push(0);
            }
            // try all possible insertions
            for l in 0..k+1 {  // l being the position to insert the job
                let mut inserted = node.inserted_jobs.clone();
                inserted.insert(l,j);
                let mut bound:Time = eqf.f[0][l] + eqf.q[0][l];
                for i in 1..m {
                    bound = max(bound, eqf.f[i as usize][l]+eqf.q[i as usize][l]);
                }
                let idletime = match self.tie_breaking {
                    TiebreakingScheme::NO => None,
                    TiebreakingScheme::FF => Some(compute_idle_time(&self.inst, &node.inserted_jobs, l, &eqf))
                };
                res.push(NEHNode {
                    inserted_jobs: inserted,
                    bound: bound,
                    idletime: idletime
                });
            }
            // compute the number of nodes with the same
            // in this code, we apply the idle time computation for each node
            // as we do not necessarily do a greedy algorithm
            // possibly, if a greedy algorithm is performed, some computations become useless
            // and the following commented code may be useful. 
            // match &self.tie_breaking {
            //     TiebreakingScheme::FF => {
            //         let mut bestbound = res[0].bound;
            //         let mut nb_best = 1;
            //         for e in res.iter() {
            //             if e.bound < bestbound {
            //                 bestbound = e.bound;
            //                 nb_best = 1;
            //             } else if e.bound == bestbound {
            //                 nb_best += 1;
            //             }
            //         }
            //         println!("nb ties: {}", nb_best);
            //     },
            //     TiebreakingScheme::NO => {}
            // }
        }
        return res;
    }
}


impl SearchTree<NEHNode, Time> for NEHSearch {
    fn root(&mut self) -> NEHNode {
        return NEHNode {
            inserted_jobs: Vec::new(),
            bound: 0,
            idletime: None,
        };
    }

    fn bound(&mut self, node: &NEHNode) -> Time {
        return node.bound;
    }

    fn goal(&mut self, node: &NEHNode) -> bool {
        return node.inserted_jobs.len() == self.inst.nb_jobs() as usize;
    }

}

impl NEHSearch {
    pub fn new(filename: &str, tb: TiebreakingScheme) -> Self {
        let inst = Instance::new(&filename).unwrap();
        // sort ordered_jobs by non-increasing sum of processing times
        let n = inst.nb_jobs();
        let m = inst.nb_machines();
        let mut ordered_jobs:Vec<JobId> = Vec::new();
        for i in 0..n {
            ordered_jobs.push(i);
        }
        let mut sum_p = vec![0;n as usize];
        for i in 0..n {
            for j in 0..m {
                sum_p[i as usize] += inst.p(i, j);
            }
        }
        ordered_jobs.sort_by_key(|a| { sum_p[*a as usize] });
        // return constructed search space
        Self {
            inst: inst,
            ordered_jobs: ordered_jobs,
            tie_breaking: tb,
        }
    }
}
