use ordered_float::OrderedFloat;

use dogs::search_space::{SearchSpace, TotalNeighborGeneration, GuidedSpace, ToSolution};

use crate::nehhelper::{compute_eqf, Eqf, compute_idle_time, LocalState};
use crate::pfsp::{JobId, Time, Instance, ProblemObjective, checker};

#[derive(Debug, Clone)]
pub enum Guide {
    Bound,  // no tie breaking
    FF,  // FF tie breaking (see https://www.sciencedirect.com/science/article/abs/pii/S0305054813003638)
    Alpha,
}


#[derive(Debug, Clone)]
pub struct NEHNode {
    /// already inserted jobs (partial sequence)
    inserted_jobs: Vec<JobId>,
    /// bound on already inserted jobs
    bound: Time,
    /// bound on remaining jobs
    bound_remaining: Time,
    // optionnal idle time of the node (used in tie breaking for instance)
    idletime: Option<Time>,
}

pub struct NEHSearch {
    inst: Instance,
    ordered_jobs: Vec<JobId>,  // ordered jobs according to the ordering heuristic
    guide: Guide,
}


impl ToSolution<NEHNode, Vec<JobId>> for NEHSearch {
    fn solution(&mut self, node: &mut NEHNode) -> Vec<JobId> {
        debug_assert!(self.goal(node));
        node.inserted_jobs.clone()
    }
}


impl SearchSpace<NEHNode, Time> for NEHSearch {
    fn initial(&mut self) -> NEHNode {
        NEHNode {
            inserted_jobs: Vec::new(),
            bound: 0,
            bound_remaining: self.inst.sum_better_first_machine()+self.inst.sum_better_last_machine(),
            idletime: None,
        }
    }

    fn bound(&mut self, node: &NEHNode) -> Time { node.bound+node.bound_remaining }

    fn g_cost(&mut self, node: &NEHNode) -> Time { node.bound }

    fn goal(&mut self, node: &NEHNode) -> bool {
        node.inserted_jobs.len() == self.inst.nb_jobs() as usize
    }

    fn handle_new_best(&mut self, mut node: NEHNode) -> NEHNode {
        // apply neighborhood decent
        let local_state = LocalState {
            s: self.solution(&mut node),
            v: self.bound(&node)
        };
        // check that the solution is correct and has the correct cost
        assert_eq!(
            checker(&self.inst, ProblemObjective::Makespan, &local_state.s),
            local_state.v
        );
        // if true {
        //     loop {
        //         match first_insertion_neighborhood(&self.inst, &local_state) {
        //             None => { break; },
        //             Some(e) => {
        //                 // println!("{} \t -> \t {}", local_state.v, e.v);
        //                 local_state = e;
        //             }
        //         }
        //     }
        // }
        // // write solution in a file
        // match &self.solution_file {
        //     None => {},
        //     Some(filename) => {
        //         let mut file = match File::create(filename) {
        //             Err(why) => panic!("couldn't create {}: {}", filename, why),
        //             Ok(file) => file
        //         };
        //         let sol = &local_state.s;
        //         for j in sol {
        //             write!(&mut file, "{} ", j).unwrap();
        //         }
        //         writeln!(&mut file, "\n").unwrap();
        //     }
        // }
        NEHNode {
            inserted_jobs: local_state.s, // already inserted jobs (partial sequence)
            bound: local_state.v,
            bound_remaining: 0,
            idletime: None,
        }
    }
}


impl GuidedSpace<NEHNode, OrderedFloat<f64>> for NEHSearch {
    fn guide(&mut self, node: &NEHNode) -> OrderedFloat<f64> {
        // let alpha:f64 = (node.inserted_jobs.len() as f64)/(self.inst.nb_jobs() as f64);
        match self.guide {
            Guide::Bound => OrderedFloat((node.bound+node.bound_remaining) as f64),
            Guide::FF => match node.idletime {
                None => OrderedFloat((node.bound+node.bound_remaining) as f64),
                Some(v) => OrderedFloat(((node.bound+node.bound_remaining) as f64)-1./(v as f64))
            }
            Guide::Alpha => {
                match node.idletime {
                    None => OrderedFloat((node.bound+node.bound_remaining) as f64),
                    Some(v) => {
                        let alpha = self.compute_alpha(node);
                        OrderedFloat(
                            alpha * (node.bound+node.bound_remaining) as f64 +
                            (1.-alpha) * v as f64
                        )
                    }
                }
            }
        }
    }
}


impl TotalNeighborGeneration<NEHNode> for NEHSearch {
    fn neighbors(&mut self, node: &mut NEHNode) -> Vec<NEHNode> {
        debug_assert!(!self.goal(node));  // should not be goal
        let mut res:Vec<NEHNode> = Vec::with_capacity(node.inserted_jobs.len()+1);
        let k = node.inserted_jobs.len();  // kth job to be inserted in the partial solution
        let j = self.ordered_jobs[k];  // job to be inserted
        let m = self.inst.nb_machines();
        // update the h-cost
        let mut new_bound_remaining = node.bound_remaining;
        if self.inst.better_on_first_machine(j) {
            new_bound_remaining -= self.inst.p(j,0);
        } else {
            new_bound_remaining -= self.inst.p(j, m-1);
        }
        if node.inserted_jobs.is_empty() {
            // if first job, just insert the first job in the ordered list
            let mut bound = 0;
            for i in 0..m {
                bound += self.inst.p(j, i);
            }
            res.push(NEHNode { 
                inserted_jobs: vec![j],
                bound,
                bound_remaining: new_bound_remaining,
                idletime: None,
            })
        } else {
            let eqf:Eqf = compute_eqf(&self.inst, &node.inserted_jobs, j);
            // try all possible insertions
            for l in 0..k+1 {  // l being the position to insert the job
                let mut inserted = node.inserted_jobs.clone();
                inserted.insert(l,j);
                let bound:Time = (0..m).map(|i|
                    eqf.f[i as usize][l] + eqf.q[i as usize][l]
                ).max().unwrap();
                let idletime = match self.guide {
                    Guide::Bound => None,
                    _ => Some(compute_idle_time(&self.inst, &node.inserted_jobs, l, &eqf)),
                };
                res.push(NEHNode {
                    inserted_jobs: inserted,
                    bound,
                    bound_remaining: new_bound_remaining,
                    idletime
                });
            }
        }
        res
    }
}


impl NEHSearch {
    pub fn new(filename: &str, guide: Guide) -> Self {
        let inst = Instance::new(filename).unwrap();
        // sort ordered_jobs by non-increasing sum of processing times
        let n = inst.nb_jobs();
        let m = inst.nb_machines();
        let mut ordered_jobs:Vec<JobId> = (0..n).collect();
        let mut sum_p:Vec<Time> = vec![0;n as usize];
        for i in 0..m {
            for j in 0..n {
                sum_p[j as usize] += inst.p(j,i);
            }
        }
        ordered_jobs.sort_by_key(|a| { -sum_p[*a as usize] });
        // return constructed search space
        Self {
            inst,
            ordered_jobs,
            guide,
        }
    }

    /**
     * computes the percentage of completion of a node:
     *  - the root is 0
     *  - a goal is 1
     *  - 0.5 indicates that there are as many jobs scheduled than unscheduled
    */
    fn compute_alpha(&self, node:&NEHNode) -> f64 {
        node.inserted_jobs.len() as f64 / self.inst.nb_jobs() as f64
    }
}
