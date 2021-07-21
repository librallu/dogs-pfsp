use bit_set::BitSet;
use ordered_float::OrderedFloat;
use std::cmp::{self, max};
use std::fs::File;
use std::io::Write;
use std::rc::Rc;
use std::cell::RefCell;

use dogs::search_space::{SearchSpace, TotalNeighborGeneration, GuidedSpace, ToSolution};
use dogs::data_structures::decision_tree::DecisionTree;
use dogs::data_structures::lazy_clonable::LazyClonable;

use crate::pfsp::{JobId, Time, Instance, MachineId, ProblemObjective, checker};
use crate::nehhelper::{LocalState, first_insertion_neighborhood};

/**
 * Guide function that helps to choose which node should be explored first
 */
#[derive(Debug, Clone)]
pub enum Guide {
    /// bound guide
    Bound,
    /// sum of idle times guide
    Idle,
    /// α * bound + (1 - α) * idle time
    /// where α: depth ratio (~0: close to the root, ~1: close to the solutions)
    Alpha,
    /// α * bound + (1 - α) * weighted idle time
    /// The weighted idle time considers more the first machine idle time on the forward part
    /// and more the last machine idle time on the backward part
    Walpha,
    /// Same as Walpha, but α is set to be the gap between the node LB and the best-known objective
    Gap,
    /// combination (by alternance of Walpha and Gap)
    Combined,
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

/**
Stores parts of the node that is costly in memory. Thus, we only use them when opening a node.
 */
#[derive(Debug, Clone)]
pub struct FBNodeLazyPart {
    /// for each machine, its first forward availability
    forward_front: Vec<Time>, 
    /// forward idle time for each machine
    forward_idle: Vec<Time>,
    /// for each machine, its first backward availability 
    backward_front: Vec<Time>,
    /// backward idle time for each machine 
    backward_idle: Vec<Time>,
    /// remaining processing time for each machine
    remaining_processing_time: Vec<Time>,
    /// subset of added jobs
    added: BitSet,
    /// total idle time of the given node
    idletime: Time,
    /// decision tree telling which choice have been done so far
    decision_tree: Rc<DecisionTree<FBDecision>>,
}

#[derive(Debug, Clone)]
pub struct FBNode {
    /// part of the node that can be computed after
    lazy_part: LazyClonable<RefCell<FBNodeLazyPart>>,
    /// bound evaluation
    bound: Time,
    /// guide evaluation
    guide: OrderedFloat<f64>,
    /// nb jobs added
    nb_added: usize,
    /// forward walpha part (idle_m / front_m)
    forward_walpha: Option<f64>,
    /// backward walpha part (idle_m / front_m)
    backward_walpha: Option<f64>,
    /// last decision taken
    last_decision: FBDecision,
}

#[derive(Debug)]
pub struct FBMakespan {
    /// instance object
    inst: Instance,
    /// best objective so far
    best_val: Option<Time>,
    /// counts number of prunings (for statistics)
    nb_prunings: u64,
    /// guide strategy used
    guide: Guide,
    /// branching scheme used
    branchingscheme: BranchingScheme,
    /// if it exists, the path to file where the solution will be contained
    solution_file: Option<String>,
    /// if true, use local search to improve found solutions 
    use_ls: bool,
    /// current iteration of the search algorithm (used to switch between guides)
    current_iter: usize,
}


impl GuidedSpace<FBNode, OrderedFloat<f64>> for FBMakespan {
    fn guide(&mut self, node: &FBNode) -> OrderedFloat<f64> { node.guide }
}


impl ToSolution<FBNode, Vec<JobId>> for FBMakespan {
    fn solution(&mut self, node: &mut FBNode) -> Vec<JobId> {
        debug_assert!(self.goal(node));
        // extracts all the decisions taken in node's branch
        let node_lazypart = node.lazy_part.lazy_get();
        let mut decisions:Vec<FBDecision> = DecisionTree::decisions_from_root(
            &node_lazypart.borrow().decision_tree
        );
        decisions.push(node.last_decision.clone());
        // build the solution from the decisions
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
}


impl SearchSpace<FBNode, Time> for FBMakespan {
    fn initial(&mut self) -> FBNode {
        let m = self.inst.nb_machines() as usize;
        // build the root node
        FBNode {
            lazy_part: LazyClonable::new(RefCell::new(FBNodeLazyPart {
                forward_front:  vec![0; m],
                backward_front: vec![0; m],
                forward_idle:   vec![0; m],
                backward_idle:  vec![0; m],
                remaining_processing_time: (0..m)
                    .map(|i| self.inst.sum_p(i as MachineId)).collect(),
                added: BitSet::new(),
                idletime: 0,
                decision_tree: DecisionTree::new(FBDecision::None),
            })),
            nb_added: 0,
            bound: self.inst.processing_time_bound(),
            guide: OrderedFloat(0.),
            forward_walpha: None,
            backward_walpha: None,
            last_decision: FBDecision::None,
        }       
    }

    fn bound(&mut self, node: &FBNode) -> Time {
        node.bound
    }

    fn g_cost(&mut self, node: &FBNode) -> Time {
        node.bound
    }

    fn goal(&mut self, node: &FBNode) -> bool { node.nb_added == self.inst.nb_jobs() as usize }

    fn handle_new_best(&mut self, mut node: FBNode) -> FBNode {
        // apply neighborhood decent
        let mut local_state = LocalState {
            s: self.solution(&mut node),
            v: self.bound(&node)
        };
        if self.use_ls {
            loop {
                match first_insertion_neighborhood(&self.inst, &local_state) {
                    None => { break; },
                    Some(e) => {
                        local_state = e;
                    }
                }
            }
        }
        // update best known solution if needed
        match self.best_val {
            None => self.best_val = Some(local_state.v),
            Some(v) => if local_state.v < v {
                self.best_val = Some(local_state.v);
            }
        }
        // check that the solution is correct and has the correct cost
        assert_eq!(
            checker(&self.inst, ProblemObjective::Makespan, &local_state.s),
            local_state.v
        );
        // write solution in a file if needed
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
        node.bound = local_state.v;
        node
    }

    fn display_statistics(&self) {
        println!();
        println!("{:>25}{:>15}", "nb pruned", self.nb_prunings);
    }

    fn restart(&mut self, _msg: String) {
        self.current_iter += 1;
    }
}


impl TotalNeighborGeneration<FBNode> for FBMakespan {
    /**
    * we assume node has its lazy part computed 
    */
    fn neighbors(&mut self, node: &mut FBNode) -> Vec<FBNode> {
        let res:Vec<FBNode>;
        match &self.branchingscheme {
            BranchingScheme::Forward => {
                res = self.generate_nodes(node, true);
            },
            BranchingScheme::Bidirectional(bidirectionaltype) => {
                match bidirectionaltype {
                    BidirectionalScheme::Alternate => {
                        if node.nb_added % 2 == 0 { // forward search
                            res = self.generate_nodes(node, true);
                        } else {  // backward search
                            res = self.generate_nodes(node, false);
                        }
                    },
                    BidirectionalScheme::MinChildren => {
                        let res_forward = self.generate_nodes(node, true);
                        let res_backward = self.generate_nodes(node, false);
                        match res_forward.len().cmp(&res_backward.len()) {
                            cmp::Ordering::Less => { res = res_forward; },
                            cmp::Ordering::Greater => { res = res_backward; },
                            cmp::Ordering::Equal => {
                                let mut sum_forward = 0;
                                for e in res_forward.iter() {
                                    sum_forward += self.bound(e);
                                }
                                let mut sum_backward = 0;
                                for e in res_backward.iter() {
                                    sum_backward += self.bound(e);
                                }
                                if sum_forward >= sum_backward {
                                    res = res_forward;
                                } else {
                                    res = res_backward;
                                }
                            },
                        }
                    }
                }
            }
        }
        res
    }
}


impl FBMakespan {
    pub fn new(filename: &str, guide:Guide, branchingscheme:BranchingScheme, use_ls:bool, solution_filename:Option<String>) -> Self {
        let inst = Instance::new(&filename).unwrap();
        Self {
            inst,
            best_val: None,
            nb_prunings: 0,
            guide,
            branchingscheme,
            solution_file: solution_filename,
            use_ls,
            current_iter: 0,
        }
    }

    /**
    use the last decision taken to update the lazy part
    */
    fn compute_lazy_part(&self, node:&mut FBNode) {
        match node.lazy_part.is_cloned() {
            true => {}  // already computed, do nothing
            false => {
                // otherwise, compute the lazy part
                let last_machine = self.inst.nb_machines()-1;
                let tmp_lazy_part = node.lazy_part.lazy_get();
                let mut res = tmp_lazy_part.borrow_mut();
                // update fronts
                match node.last_decision {
                    FBDecision::None => {},
                    FBDecision::AddForward(j) => {
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
                            // update remaining processing time
                            res.remaining_processing_time[m as usize] -= self.inst.p(j,m);            
                        }
                        // update remaining processing time for first machine
                        res.remaining_processing_time[0] -= self.inst.p(j,0);
                        res.added.insert(j as usize);
                    },
                    FBDecision::AddBackward(j) => {
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
                        // update remaining processing time for first machine
                        res.remaining_processing_time[last_machine as usize] -= self.inst.p(j,last_machine);
                        res.added.insert(j as usize);
                    }
                }
                // decision tree telling which choice have been done so far
                res.decision_tree = DecisionTree::add_child(
                    &res.decision_tree, 
                    node.last_decision.clone()
                );
            }
        }
    }

    /**
     * computes the percentage of completion of a node:
     *  - the root is 0
     *  - a goal is 1
     *  - 0.5 indicates that there are as many jobs scheduled than unscheduled
     */
    fn compute_alpha(&self, node:&FBNode) -> f64 {
        node.nb_added as f64 / self.inst.nb_jobs() as f64
    }

    fn generate_nodes(&mut self, node:&mut FBNode, generate_forward:bool) -> Vec<FBNode> {
        self.compute_lazy_part(node);  // make sure the lazy part is computed
        let mut res:Vec<FBNode> = Vec::with_capacity(self.inst.nb_jobs() as usize - node.nb_added);
        let node_lazypart = node.lazy_part.lazy_get();
        let lazy_part = node_lazypart.as_ref().borrow();
        for j in 0..self.inst.nb_jobs() {
            if !lazy_part.added.contains(j as usize) {
                let n:FBNode = self.create_child(node, j, generate_forward);
                match self.best_val {
                    None => { res.push(n); },
                    Some(v) => {
                        if self.bound(&n) >= v {
                            self.nb_prunings += 1;
                        } else {
                            res.push(n);
                        }
                    }
                }
            }
        }
        res
    }

    fn create_child(&self, node:&mut FBNode, j:JobId, generate_forward:bool) -> FBNode {
        self.compute_lazy_part(node);  // make sure the lazy part is computed
        let node_lazypart = node.lazy_part.lazy_get();
        let lazy_part = node_lazypart.as_ref().borrow();
        let mut bound = node.bound;
        let mut new_idle = 0;
        let m = self.inst.nb_machines();
        let mut new_walpha = 0.;
        // iterate over the front, update the bound
        match generate_forward {
            true => {
                let p = self.inst.p(j, 0);
                let mut current_front = lazy_part.forward_front[0] + p;
                bound = max(
                    bound, 
                    current_front + lazy_part.remaining_processing_time[0] - p + lazy_part.backward_front[0]
                );
                for i in 1..m {
                    let p = self.inst.p(j, i);
                    let start = max(current_front, lazy_part.forward_front[i as usize]);
                    let current_idle = max(0, start - lazy_part.forward_front[i as usize]);
                    new_idle += current_idle;
                    current_front = start + p;
                    bound = max(bound, current_front + lazy_part.remaining_processing_time[i as usize] - p + lazy_part.backward_front[i as usize]);
                    new_walpha += (current_idle + lazy_part.forward_idle[i as usize]) as f64 / current_front as f64;
                }
            },
            false => {
                let p = self.inst.p(j, m-1);
                let mut current_front = lazy_part.backward_front[m as usize -1] + p;
                bound = max(bound, current_front + lazy_part.remaining_processing_time[0] - p + lazy_part.forward_front[0]);
                for i in (0..m-1).rev() {
                    let p = self.inst.p(j, i);
                    let start = max(current_front, lazy_part.backward_front[i as usize]);
                    let current_idle = max(0, start - lazy_part.backward_front[i as usize]);
                    new_idle += current_idle;
                    current_front = start + p;
                    bound = max(bound, current_front + lazy_part.remaining_processing_time[i as usize] - p + lazy_part.forward_front[i as usize]);
                    new_walpha += (current_idle + lazy_part.backward_idle[i as usize]) as f64 / current_front as f64;
                }
            }
        }
        let (forward_walpha, backward_walpha) = match generate_forward {
            true  => (Some(new_walpha), node.backward_walpha),
            false => (node.forward_walpha, Some(new_walpha))
        };
        let mut res = FBNode {
            lazy_part: node.lazy_part.clone(),
            bound,
            guide: OrderedFloat(0.),
            nb_added: node.nb_added+1,
            forward_walpha,
            backward_walpha,
            last_decision: match generate_forward {
                true => FBDecision::AddForward(j),
                false => FBDecision::AddBackward(j)
            }
        };
        res.guide = match self.guide {
            Guide::Bound   => OrderedFloat(bound as f64),
            Guide::Idle    => OrderedFloat((new_idle + lazy_part.idletime) as f64),
            Guide::Alpha   => {
                let alpha = self.compute_alpha(&res);
                let n = res.nb_added as f64;
                let m = self.inst.nb_machines() as f64;
                let c = n/m;
                OrderedFloat(
                    (bound as f64) * alpha +
                    ((new_idle + lazy_part.idletime) as f64) * c * (1.-alpha)
                )
            }
            Guide::Walpha  => {
                self.create_walpha_guide(&res)
            }
            Guide::Gap => {
                self.create_gap_guide(&res)
            }
            Guide::Combined => {
                if self.current_iter % 2 == 0 {
                    self.create_walpha_guide(&res)
                } else {
                    self.create_gap_guide(&res)
                }
            }
        };
        res
    }

    fn create_walpha_guide(&self, res:&FBNode) -> OrderedFloat<f64> {
        match res.forward_walpha {
            None => OrderedFloat(1e31),
            Some(wf) => {
                match res.backward_walpha {
                    None => OrderedFloat(1e31),
                    Some(wb) => {
                        let alpha = self.compute_alpha(&res);
                        OrderedFloat(
                            alpha      * (res.bound as f64) +
                            (1.-alpha) * (wf + wb) * (res.bound as f64)
                        )
                    }
                }
            }
        }
    }

    fn create_gap_guide(&self, res:&FBNode) -> OrderedFloat<f64> {
        match self.best_val {
            None => { OrderedFloat(res.bound as f64) },
            Some(v) => {
                match res.forward_walpha {
                    None => {
                        OrderedFloat(1e31)
                    }
                    Some(wf) => {
                        match res.backward_walpha {
                            None => {
                                OrderedFloat(1e31)
                            }
                            Some(wb) => {
                                let gap:f64 = ((v-res.bound) as f64)/(v as f64);
                                // gap close to 1 : bound tight
                                // gap small (close to 0) : bound not very tight
                                OrderedFloat(
                                    (res.bound as f64)*(1./gap) + 
                                    (wf+wb)*gap
                                )
                            }
                        }
                    }
                }
            }
        }
    }
}

