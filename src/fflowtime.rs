use bit_set::BitSet;
use ordered_float::OrderedFloat;
use std::{fmt::Debug, fs::File};
use std::io::Write;
use std::rc::Rc;

use dogs::searchspace::{SearchSpace, SearchTree, TotalChildrenExpansion, GuidedSpace};
use dogs::datastructures::decisiontree::DecisionTree;

use crate::pfsp::{JobId, MachineId, Time, Instance};

pub type NodeVec = Vec<Time>;

/**
 * Guide function that helps to choose which node should be explored first
 */
#[derive(Debug, Clone)]
pub enum Guide {
    Bound,
    Idle,
    Alpha,
    Walpha,
    Walphaold
}

#[derive(Debug, Clone)]
pub enum Decision {
    /// no decision
    None,
    /// adding a job forward
    AddForwardSearch(JobId),
}

#[derive(Debug, Clone)]
pub struct Node {
    /// for each machine, its first forward availability
    forward_front: NodeVec,
    /// forward idle time for each machine
    forward_idle: NodeVec,
    /// nb jobs added
    nb_added: usize,
    /// subset of added jobs
    added: BitSet,
    /// prefix bound evaluation
    bound: Time,
    /// strong bound evaluation
    strong_bound: Time,
    /// total idle time of the given node
    idletime: Time,
    /// weighted idle time
    weightedidle: f64,
    /// decision tree telling which choice have been done so far
    decision_tree: Rc<DecisionTree<Decision>>,
}

#[derive(Debug)]
pub struct ForwardSearch {
    inst: Instance,
    guide: Guide,
    solution_file: Option<String>, // if it exists, the path to file where the solution will be contained
    decision_tree: Rc<DecisionTree<Decision>>,
}


impl GuidedSpace<Node, OrderedFloat<f64>> for ForwardSearch {
    fn guide(&mut self, node: &Node) -> OrderedFloat<f64> {
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
                let n = node.nb_added as f64;
                let m = self.inst.nb_machines() as f64;
                // let c = m.ln()/n;
                let c = n/m.ln()/n.ln();
                return OrderedFloat(
                    alpha * (self.bound(node) as f64) +
                    (1.-alpha) * c * node.weightedidle
                );
            },
            Guide::Walphaold => {
                let alpha = self.compute_alpha(node);
                let bound_part:f64 = self.bound(node) as f64;
                let widle_part:f64 = (node.weightedidle as f64 + self.inst.nb_machines() as f64 * node.idletime as f64) as f64 / 2.;
                return OrderedFloat(alpha * bound_part + (1.-alpha) * widle_part);
            }
        }
    }
}


impl SearchSpace<Node, Vec<JobId>> for ForwardSearch {
    fn solution(&mut self, node: &Node) -> Vec<JobId> {
        debug_assert!(self.goal(node));
        let decisions:Vec<Decision> = DecisionTree::decisions_from_root(&node.decision_tree);
        let mut res:Vec<JobId> = Vec::with_capacity(self.inst.nb_jobs() as usize);
        for e in decisions {
            match e {
                Decision::AddForwardSearch(j) => { res.push(j); }
                _ => {}
            }
        }
        res
    }

    fn handle_new_best(&mut self, node: Node) -> Node {
        let sol = self.solution(&node);
        // write solution in a file
        match &self.solution_file {
            None => {},
            Some(filename) => {
                let mut file = match File::create(filename) {
                    Err(why) => panic!("couldn't create {}: {}", filename, why),
                    Ok(file) => file
                };
                for j in sol {
                    write!(&mut file, "{} ", j).unwrap();
                }
                writeln!(&mut file, "\n").unwrap();
            }
        }
        node
    }
}


impl TotalChildrenExpansion<Node> for ForwardSearch {
    fn children(&mut self, node: &mut Node) -> Vec<Node> {
        let mut res = Vec::with_capacity(self.inst.nb_jobs() as usize - node.nb_added);
        for j in 0..self.inst.nb_jobs() {
            if !node.added.contains(j as usize) {
                res.push(self.add_job_forward(node, j));
            }
        }
        res
    }
}


impl SearchTree<Node, Time> for ForwardSearch {

    fn root(&mut self) -> Node {
        let m = self.inst.nb_machines() as usize;
        Node {
            forward_front: vec![0; m],
            forward_idle: vec![0; m],
            nb_added: 0,
            added: BitSet::new(),
            strong_bound: self.inst.sum_p((m-1) as MachineId),
            bound: 0,
            idletime: 0,
            weightedidle: 0.,
            decision_tree: self.decision_tree.clone(),
        }
    }

    fn bound(&mut self, node: &Node) -> Time {
        return node.bound;
    }

    fn goal(&mut self, node: &Node) -> bool {
        return node.nb_added == self.inst.nb_jobs() as usize;
    }
}

impl ForwardSearch {
    pub fn new(filename: &str, guide:Guide, solution_filename:Option<String>) -> Self {
        Self {
            inst: Instance::new(&filename).unwrap(),
            guide,
            solution_file: solution_filename,
            decision_tree: DecisionTree::new(Decision::None),
        }
    }
    

    pub fn add_job_forward(&self, node:&Node, j: JobId) -> Node {
        let mut res = node.clone();
        let alpha = self.compute_alpha(node)+1.;
        // update front & idle time
        // let n = node.nb_added as f64 + 1.;
        res.forward_front[0] += self.inst.p(j, 0);
        for m in 1..self.inst.nb_machines() {
            let start:Time;
            if res.forward_front[(m-1) as usize] > res.forward_front[m as usize] {
                start = res.forward_front[(m-1) as usize];
                let new_idle_time = start - res.forward_front[m as usize];
                res.forward_idle[m as usize] += new_idle_time;
                res.idletime += new_idle_time;
                res.weightedidle += new_idle_time as f64 * (alpha*(self.inst.nb_machines()-m) as f64);
                // res.weightedidle += new_idle_time as f64 * (1./m as f64);
            } else {
                start = res.forward_front[m as usize];
            }
            res.forward_front[m as usize] = start + self.inst.p(j,m);
        }
        // register that the job is added
        res.decision_tree = DecisionTree::add_child(&node.decision_tree, Decision::AddForwardSearch(j));
        res.added.insert(j as usize);
        res.nb_added += 1;
        // update the bound
        let lastmachineid = (self.inst.nb_machines()-1) as usize;
        res.bound += res.forward_front[lastmachineid];
        return res;
    }

    /**
     * computes the percentage of completion of a node:
     *  - the root is 0
     *  - a goal is 1
     *  - 0.5 indicates that there are as many jobs scheduled than unscheduled
     */
    fn compute_alpha(&self, node:&Node) -> f64 {
        ( node.nb_added as f64 ) / ( self.inst.nb_jobs() as f64 )
    }
}

