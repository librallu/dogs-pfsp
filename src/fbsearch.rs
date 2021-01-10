use bit_set::BitSet;
use ordered_float::OrderedFloat;
use std::cmp::max;
use std::fs::File;
use std::io::Write;

use dogs::searchspace::{SearchSpace, SearchTree, TotalChildrenExpansion, GuidedSpace};

use crate::pfsp::{JobId, Time, Instance};

/**
 * Objective function to minimize
 */
#[derive(Debug)]
pub enum Objective {
    Makespan,
    Flowtime
}

/**
 * Guide function that helps to choose which node should be explored first
 */
#[derive(Debug)]
pub enum Guide {
    Bound, // bound
    IdleTime, // idle
    BoundAlphaIdle, // alpha
    BoundWAlphaIdle,  // wfrontalpha
    BoundW2AlphaIdle,  // walpha
    BoundGapIdle, // gap
}

/**
 * type of bidirectional search used
 */
#[derive(Debug)]
pub enum BidirectionalScheme {
    Alternate, // alternates between forward and backward
    MinChildren, // chooses the direction which has the less children (larger sum of bounds if ties happen)
}

/**
 * Ways to expand a node
 */
#[derive(Debug)]
pub enum BranchingScheme {
    Forward,
    Bidirectional(BidirectionalScheme), 
}

#[derive(Debug, Clone)]
pub struct FBNode {
    forward_jobs: Vec<JobId>,  // partial job order
    backward_jobs: Vec<JobId>,  // partial job order
    forward_front: Vec<Time>,  // for each machine, its first forward availability
    backward_front: Vec<Time>,  // for each machine, its first backward availability
    forward_idle: Vec<Time>,  // forward idle time for each machine
    backward_idle: Vec<Time>, // backward idle time for each machine 
    remaining_processing_time: Vec<Time>,  // remaining processing time for each machine
    added: BitSet,  // subset of added jobs
    bound: Time,  // bound evaluation
    idletime: Time,  // total idle time of the given node
    weightedidle: f64, // weighted idle time
}

#[derive(Debug)]
pub struct FBSearch {
    inst: Instance,
    best_val: Option<Time>,
    nb_prunings: u64,
    objtype: Objective,
    guide: Guide,
    branchingscheme: BranchingScheme,
    solution_file: Option<String>, // if it exists, the path to file where the solution will be contained 
}


impl GuidedSpace<FBNode, OrderedFloat<f64>> for FBSearch {
    fn guide(&mut self, node: &FBNode) -> OrderedFloat<f64> {
        match self.guide {
            Guide::Bound => { return OrderedFloat(node.bound as f64); },
            Guide::IdleTime => { return OrderedFloat(node.idletime as f64); },
            Guide::BoundAlphaIdle => {
                let alpha = self.compute_alpha(node);
                let c = (self.inst.nb_jobs()/self.inst.nb_machines()) as f64;
                return OrderedFloat(
                    alpha*(node.bound as f64) +
                    (1.-alpha)*c*(node.idletime as f64)
                );
            },
            Guide::BoundWAlphaIdle => {
                let alpha = self.compute_alpha(node);
                match self.compute_wfront(&node) {
                    None => { return OrderedFloat(1e31); }
                    Some(widle) => {
                        return OrderedFloat(
                            alpha*(node.bound as f64) +
                            (1.-alpha)*widle*(node.bound as f64)
                        );
                    }
                }
            },
            Guide::BoundGapIdle => {
                match self.best_val {
                    None => { return OrderedFloat(node.bound as f64); },
                    Some(v) => {
                        match self.compute_wfront(&node) {
                            None => { return OrderedFloat(1e31); } // infinity
                            Some(w_front_idle) => {
                                let gap:f64 = ((v-node.bound) as f64)/(v as f64);
                                // gap close to 1 : bound tight
                                // gap small (close to 0) : bound not very tight
                                return OrderedFloat(
                                    (node.bound as f64)*(1./gap)+
                                    w_front_idle*gap
                                );
                            }
                        }
                    }
                }
            },
            Guide::BoundW2AlphaIdle => {
                let alpha = self.compute_alpha(node);
                let bound_part:f64 = self.bound(node) as f64;
                let widle_part:f64 = (node.weightedidle as f64 + self.inst.nb_machines() as f64 * node.idletime as f64) as f64 / 2.;
                return OrderedFloat(alpha * bound_part + (1.-alpha) * widle_part);
            }
        }
    }
}


impl SearchSpace<FBNode, Vec<JobId>> for FBSearch {
    fn solution(&mut self, node: &FBNode) -> Vec<JobId> {
        debug_assert!(self.goal(node));
        let mut res:Vec<JobId> = Vec::new();
        for e in node.forward_jobs.iter() {
            res.push(*e);
        }
        for e in node.backward_jobs.iter().rev() {
            res.push(*e);
        }
        return res;
    }

    fn handle_new_best(&mut self, node: &FBNode) {
        let new_val = self.bound(node);
        match self.best_val {
            None => self.best_val = Some(new_val),
            Some(v) => if new_val < v {
                self.best_val = Some(new_val);
            }
        }
        // write solution in a file
        match &self.solution_file {
            None => {},
            Some(filename) => {
                let mut file = match File::create(filename) {
                    Err(why) => panic!("couldn't create {}: {}", filename, why),
                    Ok(file) => file
                };
                let sol = self.solution(node);
                for j in sol {
                    write!(&mut file, "{} ", j).unwrap();
                }
                writeln!(&mut file, "\n").unwrap();
            }
        }
    }

    fn display_statistics(&self) {
        println!();
        println!("{:>25}{:>15}", "nb pruned", self.nb_prunings);
    }
}


impl TotalChildrenExpansion<FBNode> for FBSearch {
    fn children(&mut self, node: &mut FBNode) -> Vec<FBNode> {
        let res:Vec<FBNode>;
        match &self.branchingscheme {
            BranchingScheme::Forward => {
                res = self.generate_forward_nodes(node);
            },
            BranchingScheme::Bidirectional(bidirectionaltype) => {
                match bidirectionaltype {
                    BidirectionalScheme::Alternate => {
                        if node.forward_jobs.len() + node.backward_jobs.len() % 2 == 0 { // forward search
                            res = self.generate_forward_nodes(node);
                        } else {  // backward search
                            res = self.generate_backward_nodes(node);
                        }
                    },
                    BidirectionalScheme::MinChildren => {
                        let res_forward = self.generate_forward_nodes(node);
                        let res_backward = self.generate_backward_nodes(node);
                        if res_forward.len() < res_backward.len() { // do forward search
                            res = res_forward;
                        } else if res_backward.len() < res_forward.len() {  // do backward search
                            res = res_backward;
                        } else {  // break ties (compute sum of bounds and chose the largest)
                            let mut sum_forward = 0;
                            for e in res_forward.iter() {
                                sum_forward += e.bound;
                            }
                            let mut sum_backward = 0;
                            for e in res_backward.iter() {
                                sum_backward += e.bound;
                            }
                            if sum_forward >= sum_backward {
                                res = res_forward;
                            } else {
                                res = res_backward;
                            }
                        }
                    }
                }
            }
        }
        return res;
    }
}


impl SearchTree<FBNode, Time> for FBSearch {

    fn root(&mut self) -> FBNode {
        let m = self.inst.nb_machines() as usize;
        let mut sum_p = vec![0; m];
        let mut bound = 0;
        match &self.objtype {
            Objective::Makespan => {
                for i in 0..self.inst.nb_machines() {
                    for j in 0..self.inst.nb_jobs() {
                        sum_p[i as usize] += self.inst.p(j,i);
                    }
                    bound = max(bound, sum_p[i as usize]);
                }
            },
            Objective::Flowtime => {}  // bound is already to zero
        }

        return FBNode {
            forward_jobs: Vec::new(),
            backward_jobs: Vec::new(),
            forward_front: vec![0; m],
            backward_front: vec![0; m],
            forward_idle: vec![0; m],
            backward_idle: vec![0; m],
            remaining_processing_time: sum_p,
            added: BitSet::new(),
            bound: bound,
            idletime: 0,
            weightedidle: 0.,
        };
    }

    fn bound(&mut self, node: &FBNode) -> Time {
        return node.bound;
    }

    fn goal(&mut self, node: &FBNode) -> bool {
        return (node.forward_jobs.len() + node.backward_jobs.len()) == self.inst.nb_jobs() as usize;
    }
}

impl FBSearch {
    pub fn new(filename: &str, objtype:Objective, guide:Guide, branchingscheme:BranchingScheme, solution_filename:Option<String>) -> Self {
        Self {
            inst: Instance::new(&filename).unwrap(),
            best_val: None,
            nb_prunings: 0,
            objtype: objtype,
            guide: guide,
            branchingscheme: branchingscheme,
            solution_file: solution_filename,
        }
    }
    

    pub fn add_job_forward(&self, node:&FBNode, j: JobId) -> FBNode {
        let mut res = node.clone();
        let last_machine = self.inst.nb_machines()-1;
        let alpha = self.compute_alpha(node)+1.;
        // update front
        res.forward_front[0] += self.inst.p(j, 0);
        for m in 1..self.inst.nb_machines() {
            let start:Time;
            if res.forward_front[(m-1) as usize] > res.forward_front[m as usize] {
                start = res.forward_front[(m-1) as usize];
                let new_idle_time = start - res.forward_front[m as usize];
                res.forward_idle[m as usize] += new_idle_time;
                res.idletime += new_idle_time;
                res.weightedidle += new_idle_time as f64 * (alpha*(self.inst.nb_machines()-m) as f64);
            } else {
                start = res.forward_front[m as usize];
            }
            res.forward_front[m as usize] = start + self.inst.p(j,m);
            res.remaining_processing_time[m as usize] -= self.inst.p(j,m); // update remaining processing time
        }
        res.remaining_processing_time[0] -= self.inst.p(j,0); // update remaining processing time for first machine
        // register that the job is added
        res.forward_jobs.push(j);
        res.added.insert(j as usize);
        // update the bound
        match self.objtype {
            Objective::Makespan => {
                match &self.branchingscheme {
                    BranchingScheme::Forward => {
                        res.bound = res.forward_front[last_machine as usize] + res.remaining_processing_time[last_machine as usize];
                    },
                    BranchingScheme::Bidirectional(_) => {
                        for m in 0..self.inst.nb_machines() {
                            res.bound = max(res.bound, 
                                res.forward_front[m as usize]+
                                res.remaining_processing_time[m as usize]+
                                res.backward_front[m as usize]
                            );
                        }
                    }
                }
            },
            Objective::Flowtime => {
                match self.branchingscheme {
                    BranchingScheme::Forward => {
                        let lastmachineid = (self.inst.nb_machines()-1) as usize;
                        res.bound += res.forward_front[lastmachineid];
                    },
                    BranchingScheme::Bidirectional(_) => {
                        panic!("bidirectional search not implemented for the flowtime minimization.");
                    }
                };
            }
        }
        return res;
    }

    pub fn add_job_backward(&self, node:&FBNode, j: JobId) -> FBNode {
        let mut res = node.clone();
        let last_machine = self.inst.nb_machines()-1;
        let alpha = self.compute_alpha(node);
        // update front
        res.backward_front[last_machine as usize] += self.inst.p(j, last_machine);
        for m in (0..last_machine).rev() {
            let start:Time;
            if res.backward_front[(m+1) as usize] > res.backward_front[m as usize] {
                start = res.backward_front[(m+1) as usize];
                let new_idle_time = res.backward_front[(m+1) as usize] - res.backward_front[m as usize];
                res.backward_idle[m as usize] += new_idle_time;
                res.idletime += new_idle_time;
                res.weightedidle += new_idle_time as f64 * (alpha*m as f64+1.);
            } else {
                start = res.backward_front[m as usize];
            }
            res.backward_front[m as usize] = start + self.inst.p(j,m);
            res.remaining_processing_time[m as usize] -= self.inst.p(j,m); // update remaining processing time
        }
        res.remaining_processing_time[last_machine as usize] -= self.inst.p(j,last_machine); // update remaining processing time for first machine
        // register node added
        res.backward_jobs.push(j);
        res.added.insert(j as usize);
        // update bound
        match self.objtype {
            Objective::Flowtime => { panic!("bidirectional search not implemented for the flowtime minimization."); },
            Objective::Makespan => {
                match self.branchingscheme {
                    BranchingScheme::Forward => { panic!("adding backwards should not be done in a forward search."); },
                    BranchingScheme::Bidirectional(_) => {
                        for m in 0..self.inst.nb_machines() {
                            res.bound = max(res.bound,
                                res.forward_front[m as usize]+
                                res.remaining_processing_time[m as usize]+
                                res.backward_front[m as usize]
                            );
                        }
                    }
                }
            }
        }
        return res;
    }

    /**
     * computes the percentage of completion of a node:
     *  - the root is 0
     *  - a goal is 1
     *  - 0.5 indicates that there are as many jobs scheduled than unscheduled
     */
    fn compute_alpha(&self, node:&FBNode) -> f64 {
        return ( node.forward_jobs.len() as f64 + node.backward_jobs.len() as f64 ) / ( self.inst.nb_jobs() as f64 );
    }

    fn generate_forward_nodes(&mut self, node:&FBNode) -> Vec<FBNode> {
        let mut res = Vec::new();
        for j in 0..self.inst.nb_jobs() {
            if !node.added.contains(j as usize) {
                let n = self.add_job_forward(node, j);
                match self.best_val {
                    None => { res.push(n); },
                    Some(v) => {
                        if n.bound >= v {
                            self.nb_prunings += 1;
                        } else {
                            res.push(n);
                        }
                    }
                }
            }
        }
        return res;
    }

    fn generate_backward_nodes(&mut self, node:&FBNode) -> Vec<FBNode> {
        let mut res = Vec::new();
        for j in 0..self.inst.nb_jobs() {
            if !node.added.contains(j as usize) {
                let n = self.add_job_backward(node, j);
                match self.best_val {
                    None => { res.push(n); },
                    Some(v) => {
                        if n.bound >= v {
                            self.nb_prunings += 1;
                        } else {
                            res.push(n);
                        }
                    }
                }
            }
        }
        return res;
    }

    fn compute_wfront(&self, node:&FBNode) -> Option<f64> {
        let mut wsum_idle = 0.;
        for i in 1..self.inst.nb_machines() {
            if node.forward_front[i as usize] == 0 {
                return None;
            }
            wsum_idle += (node.forward_idle[i as usize] as f64) / (node.forward_front[i as usize] as f64);
        }
        for i in 0..self.inst.nb_machines()-1 {
            if node.backward_front[i as usize] == 0 {
                return None;
            }
            wsum_idle += (node.backward_idle[i as usize] as f64) / (node.backward_front[i as usize] as f64);
        }
        return Some(wsum_idle);
    }
}

