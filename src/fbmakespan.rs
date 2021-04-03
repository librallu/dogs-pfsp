use bit_set::BitSet;
use ordered_float::OrderedFloat;
use std::cmp::max;
use std::fs::File;
use std::io::Write;
use std::rc::Rc;

use dogs::searchspace::{SearchSpace, SearchTree, TotalChildrenExpansion, GuidedSpace};
use dogs::datastructures::decisiontree::DecisionTree;

use crate::pfsp::{JobId, MachineId, Time, Instance};
use crate::nehhelper::{LocalState, first_insertion_neighborhood};

/**
 * Guide function that helps to choose which node should be explored first
 */
#[derive(Debug, Clone)]
pub enum Guide {
    Bound, // bound
    Idle, // idle
    Alpha, // alpha
    Walpha,  // wfrontalpha
    Gap, // gap
}

/**
 * type of bidirectional search used
 */
#[derive(Debug, Clone)]
pub enum BidirectionalScheme {
    /// alternates between forward and backward
    Alternate,
    /// chooses the direction which has the less children (larger sum of bounds if ties happen)
    MinChildren,
}

#[derive(Debug, Clone)]
pub enum FBDecision {
    /// no decision
    None,
    /// adding a job forward
    AddForward(JobId),
    /// adding a job backward
    AddBackward(JobId)
}

/**
 * Ways to expand a node
 */
#[derive(Debug, Clone)]
pub enum BranchingScheme {
    Forward,
    Bidirectional(BidirectionalScheme), 
}

#[derive(Debug, Clone)]
pub struct FBNode {
    forward_front: Vec<Time>,  // for each machine, its first forward availability
    backward_front: Vec<Time>,  // for each machine, its first backward availability
    forward_idle: Vec<Time>,  // forward idle time for each machine
    backward_idle: Vec<Time>, // backward idle time for each machine 
    remaining_processing_time: Vec<Time>,  // remaining processing time for each machine
    nb_added: usize, /// nb jobs added
    added: BitSet,  // subset of added jobs
    bound: Time,  // bound evaluation
    idletime: Time,  // total idle time of the given node
    /// decision tree telling which choice have been done so far
    decision_tree: Rc<DecisionTree<FBDecision>>,
}

#[derive(Debug)]
pub struct FBMakespan {
    inst: Instance,
    best_val: Option<Time>,
    nb_prunings: u64,
    guide: Guide,
    branchingscheme: BranchingScheme,
    solution_file: Option<String>, // if it exists, the path to file where the solution will be contained
    use_ls: bool, // if true, use local search to improve found solutions
    decision_tree: Rc<DecisionTree<FBDecision>>,
}


impl GuidedSpace<FBNode, OrderedFloat<f64>> for FBMakespan {
    fn guide(&mut self, node: &FBNode) -> OrderedFloat<f64> {
        match self.guide {
            Guide::Bound => { return OrderedFloat(node.bound as f64); },
            Guide::Idle => { return OrderedFloat(node.idletime as f64); },
            Guide::Alpha => {
                let alpha = self.compute_alpha(node);
                let n = node.nb_added as f64;
                let m = self.inst.nb_machines() as f64;
                let c = n/m;
                return OrderedFloat(
                    alpha*(node.bound as f64) +
                    (1.-alpha)*c*(node.idletime as f64)
                );
            },
            Guide::Walpha => {
                let alpha = self.compute_alpha(node);
                match self.compute_wfront(&node) {
                    None => { return OrderedFloat(1e31); }
                    Some(widle) => {
                        return OrderedFloat(
                            alpha * (node.bound as f64) +
                            (1.-alpha) * widle * (node.bound as f64)
                        );
                    }
                }
            },
            Guide::Gap => {
                match self.best_val {
                    None => { return OrderedFloat(node.bound as f64); },
                    Some(v) => {
                        match self.compute_wfront(&node) {
                            None => { return OrderedFloat(1e31); } // infinity because only one direction
                            Some(w_front_idle) => {
                                let gap:f64 = ((v-node.bound) as f64)/(v as f64);
                                // gap close to 1 : bound tight
                                // gap small (close to 0) : bound not very tight
                                return OrderedFloat(
                                    (node.bound as f64)*(1./gap) + 
                                    w_front_idle*gap
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}


impl SearchSpace<FBNode, Vec<JobId>> for FBMakespan {
    fn solution(&mut self, node: &FBNode) -> Vec<JobId> {
        debug_assert!(self.goal(node));
        let decisions:Vec<FBDecision> = DecisionTree::decisions_from_root(&node.decision_tree);
        let mut res:Vec<JobId> = Vec::with_capacity(self.inst.nb_jobs() as usize);
        let mut jobs_backward:Vec<JobId> = Vec::new();
        for e in decisions {
            match e {
                FBDecision::AddForward(j) => { res.push(j); }
                FBDecision::AddBackward(j) => { jobs_backward.push(j); }
                _ => {}
            }
        }
        jobs_backward.reverse();
        res.append(&mut jobs_backward);
        res
    }

    fn handle_new_best(&mut self, node: FBNode) -> FBNode {
        // apply neighborhood decent
        let mut local_state = LocalState {
            s: self.solution(&node),
            v: self.bound(&node)
        };
        if self.use_ls {
            loop {
                match first_insertion_neighborhood(&self.inst, &local_state) {
                    None => { break; },
                    Some(e) => {
                        // println!("{} \t -> \t {}", local_state.v, e.v);
                        local_state = e;
                    }
                }
            }
        }
        // update best known solution
        match self.best_val {
            None => self.best_val = Some(local_state.v),
            Some(v) => if local_state.v < v {
                self.best_val = Some(local_state.v);
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
                let sol = &local_state.s;
                for j in sol {
                    write!(&mut file, "{} ", j).unwrap();
                }
                writeln!(&mut file, "\n").unwrap();
            }
        }
        let mut res = node.clone();
        res.bound = local_state.v;
        res
    }

    fn display_statistics(&self) {
        println!();
        println!("{:>25}{:>15}", "nb pruned", self.nb_prunings);
    }
}


impl TotalChildrenExpansion<FBNode> for FBMakespan {
    fn children(&mut self, node: &mut FBNode) -> Vec<FBNode> {
        match &self.branchingscheme {
            BranchingScheme::Forward => {
                self.generate_forward_nodes(node)
            },
            BranchingScheme::Bidirectional(bidirectionaltype) => {
                match bidirectionaltype {
                    BidirectionalScheme::Alternate => {
                        if node.nb_added % 2 == 0 { // forward search
                            self.generate_forward_nodes(node)
                        } else {  // backward search
                            self.generate_backward_nodes(node)
                        }
                    },
                    BidirectionalScheme::MinChildren => {
                        let res_forward = self.generate_forward_nodes(node);
                        let res_backward = self.generate_backward_nodes(node);
                        if res_forward.len() < res_backward.len() { // do forward search
                            res_forward
                        } else if res_backward.len() < res_forward.len() {  // do backward search
                            res_backward
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
                                res_forward
                            } else {
                                res_backward
                            }
                        }
                    }
                }
            }
        }
    }
}


impl SearchTree<FBNode, Time> for FBMakespan {

    fn root(&mut self) -> FBNode {
        let m = self.inst.nb_machines() as usize;
        let mut sum_p ;
        let mut bound = 0;
        sum_p = vec![0; m];
        for i in 0..self.inst.nb_machines() {
            sum_p[i as usize] = self.inst.sum_p(i as MachineId);
            bound = max(bound, sum_p[i as usize]);
        }
        let (backward_front,backward_idle) = match self.branchingscheme {
            BranchingScheme::Forward => { (Vec::new(), Vec::new()) }
            BranchingScheme::Bidirectional(_) => { (vec![0; m],vec![0; m]) }
        };
        FBNode {
            forward_front: vec![0; m],
            backward_front,
            forward_idle: vec![0; m],
            backward_idle,
            remaining_processing_time: sum_p,
            nb_added: 0,
            added: BitSet::new(),
            bound,
            idletime: 0,
            decision_tree: self.decision_tree.clone(),
        }
    }

    fn bound(&mut self, node: &FBNode) -> Time {
        return node.bound;
    }

    fn goal(&mut self, node: &FBNode) -> bool {
        return node.nb_added == self.inst.nb_jobs() as usize;
    }
}

impl FBMakespan {
    pub fn new(filename: &str, guide:Guide, branchingscheme:BranchingScheme, use_ls:bool, solution_filename:Option<String>) -> Self {
        Self {
            inst: Instance::new(&filename).unwrap(),
            best_val: None,
            nb_prunings: 0,
            guide,
            branchingscheme,
            solution_file: solution_filename,
            use_ls,
            decision_tree: DecisionTree::new(FBDecision::None),
        }
    }
    

    pub fn add_job_forward(&self, node:&FBNode, j: JobId) -> FBNode {
        let mut res = node.clone();
        let last_machine = self.inst.nb_machines()-1;
        // update front
        res.forward_front[0] += self.inst.p(j, 0);
        for m in 1..self.inst.nb_machines() {
            let start:Time;
            if res.forward_front[(m-1) as usize] > res.forward_front[m as usize] {
                start = res.forward_front[(m-1) as usize];
                let new_idle_time = start - res.forward_front[m as usize];
                res.forward_idle[m as usize] += new_idle_time;
                res.idletime += new_idle_time;
            } else {
                start = res.forward_front[m as usize];
            }
            res.forward_front[m as usize] = start + self.inst.p(j,m);
            // update remaining processing time (only for makespan)
            res.remaining_processing_time[m as usize] -= self.inst.p(j,m);            
        }
        // update remaining processing time for first machine
        res.remaining_processing_time[0] -= self.inst.p(j,0);
        // register that the job is added
        res.decision_tree = DecisionTree::add_child(&node.decision_tree, FBDecision::AddForward(j));
        res.added.insert(j as usize);
        res.nb_added += 1;
        // update the bound
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
        return res;
    }

    pub fn add_job_backward(&self, node:&FBNode, j: JobId) -> FBNode {
        let mut res = node.clone();
        let last_machine = self.inst.nb_machines()-1;
        // update front
        res.backward_front[last_machine as usize] += self.inst.p(j, last_machine);
        for m in (0..last_machine).rev() {
            let start:Time;
            if res.backward_front[(m+1) as usize] > res.backward_front[m as usize] {
                start = res.backward_front[(m+1) as usize];
                let new_idle_time = res.backward_front[(m+1) as usize] - res.backward_front[m as usize];
                res.backward_idle[m as usize] += new_idle_time;
                res.idletime += new_idle_time;
            } else {
                start = res.backward_front[m as usize];
            }
            res.backward_front[m as usize] = start + self.inst.p(j,m);
            res.remaining_processing_time[m as usize] -= self.inst.p(j,m); // update remaining processing time
        }
        res.remaining_processing_time[last_machine as usize] -= self.inst.p(j,last_machine); // update remaining processing time for first machine
        // register node added
        res.decision_tree = DecisionTree::add_child(&node.decision_tree, FBDecision::AddBackward(j));
        res.added.insert(j as usize);
        res.nb_added += 1;
        // update bound
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
        return res;
    }

    /**
     * computes the percentage of completion of a node:
     *  - the root is 0
     *  - a goal is 1
     *  - 0.5 indicates that there are as many jobs scheduled than unscheduled
     */
    fn compute_alpha(&self, node:&FBNode) -> f64 {
        return ( node.nb_added as f64 ) / ( self.inst.nb_jobs() as f64 );
    }

    fn generate_forward_nodes(&mut self, node:&FBNode) -> Vec<FBNode> {
        let mut res = Vec::with_capacity(self.inst.nb_jobs() as usize - node.nb_added);
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
        let mut res = Vec::with_capacity(self.inst.nb_jobs() as usize - node.nb_added);
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

