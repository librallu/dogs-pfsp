use bit_set::BitSet;
use ordered_float::OrderedFloat;
use std::{fmt::Debug, fs::File};
use std::io::Write;
use std::rc::Rc;
use std::cell::RefCell;
use std::cmp::max;

use dogs::search_space::{SearchSpace, TotalNeighborGeneration, GuidedSpace, ToSolution};
use dogs::data_structures::decision_tree::DecisionTree;
use dogs::data_structures::lazy_clonable::LazyClonable;

use crate::pfsp::{JobId, Time, Instance, ProblemObjective, checker};

pub type NodeVec = Vec<Time>;

/**
 * Guide function that helps to choose which node should be explored first
 */
#[derive(Debug, Clone)]
pub enum Guide {
    Bound,
    Idle,
    Alpha,
}

#[derive(Debug, Clone)]
pub enum Decision {
    /// no decision
    None,
    /// adding a job forward
    AddForwardSearch(JobId),
}

#[derive(Debug, Clone)]
pub struct NodeLazyPart {
    /// for each machine, its first forward availability
    forward_front: NodeVec,
    /// subset of added jobs
    added: BitSet,
    /// decision tree telling which choice have been done so far
    decision_tree: Rc<DecisionTree<Decision>>,
}

#[derive(Debug, Clone)]
pub struct Node {
    /// nb jobs added
    nb_added: usize,
    /// prefix bound evaluation (g-cost)
    bound: Time,
    /// h-cost estimate
    fcost: Time,
    /// total idle time of the given node
    idletime: Time,
    /// weighted idle time
    weightedidle: f64,
    /// last decision
    decision: Decision,
    /// lazy part of the node
    lazy_part: LazyClonable<RefCell<NodeLazyPart>>,
}

#[derive(Debug)]
pub struct ForwardSearch {
    inst: Instance,
    guide: Guide,
    solution_file: Option<String>, // if it exists, the path to file where the solution will be contained
    decision_tree: Rc<DecisionTree<Decision>>,}


impl GuidedSpace<Node, OrderedFloat<f64>> for ForwardSearch {
    fn guide(&mut self, node: &Node) -> OrderedFloat<f64> {
        match self.guide {
            Guide::Bound => { OrderedFloat(node.bound as f64) },
            Guide::Idle => { OrderedFloat(node.idletime as f64) },
            Guide::Alpha => {
                let alpha = self.compute_alpha(node);
                let c = node.nb_added as f64 / self.inst.nb_machines() as f64;
                OrderedFloat(
                    alpha*(node.bound as f64) +
                    (1.-alpha)*c*(node.idletime as f64)
                )
            },
        }
    }
}


impl ToSolution<Node, Vec<JobId>> for ForwardSearch {
    fn solution(&mut self, node: &mut Node) -> Vec<JobId> {
        debug_assert!(self.goal(node));
        self.compute_lazy_part(node);  // make sure the lazy part is computed
        let node_lazypart = node.lazy_part.lazy_get();
        let lazy_part = node_lazypart.as_ref().borrow();
        let decisions:Vec<Decision> = DecisionTree::decisions_from_root(&lazy_part.decision_tree);
        let mut res:Vec<JobId> = Vec::with_capacity(self.inst.nb_jobs() as usize);
        for e in decisions {
            if let Decision::AddForwardSearch(j) = e { res.push(j); }
        }
        res
    }
}


impl SearchSpace<Node, Time> for ForwardSearch {

    fn initial(&mut self) -> Node {
        let m = self.inst.nb_machines() as usize;
        Node {
            nb_added: 0,
            bound: 0,
            fcost: 0,
            idletime: 0,
            weightedidle: 0.,
            decision: Decision::None,
            lazy_part: LazyClonable::new(RefCell::new(NodeLazyPart {
                forward_front: vec![0; m],
                added: BitSet::new(),
                decision_tree: self.decision_tree.clone(),
            })),
        }
    }

    fn g_cost(&mut self, node: &Node) -> Time { node.bound }

    fn bound(&mut self, node: &Node) -> Time { node.fcost }

    fn goal(&mut self, node: &Node) -> bool { node.nb_added == self.inst.nb_jobs() as usize }

    fn handle_new_best(&mut self, mut node: Node) -> Node {
        let sol = self.solution(&mut node);
        // check that the solution is correct and has the correct cost
        assert_eq!(
            checker(&self.inst, ProblemObjective::Flowtime, &sol),
            node.bound
        );
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


impl TotalNeighborGeneration<Node> for ForwardSearch {
    fn neighbors(&mut self, node: &mut Node) -> Vec<Node> {
        self.compute_lazy_part(node);  // make sure the lazy part is computed
        let node_lazypart = node.lazy_part.lazy_get();
        let lazy_part = node_lazypart.as_ref().borrow();
        let mut res = Vec::with_capacity(self.inst.nb_jobs() as usize - node.nb_added);
        for j in 0..self.inst.nb_jobs() {
            if !lazy_part.added.contains(j as usize) {
                res.push(self.add_job_forward(node, j));
            }
        }
        res
    }
}


impl ForwardSearch {
    pub fn new(filename: &str, guide:Guide, solution_filename:Option<String>) -> Self {
        let inst = Instance::new(&filename).unwrap();
        Self {
            inst,
            guide,
            solution_file: solution_filename,
            decision_tree: DecisionTree::new(Decision::None),
        }
    }


    pub fn compute_lazy_part(&self, node:&mut Node) {
        match node.lazy_part.is_cloned() {
            true => {}  // already computed, do nothing
            false => {
                let node_lazypart:Rc<RefCell<NodeLazyPart>> = node.lazy_part.lazy_get();
                let mut res = node_lazypart.as_ref().borrow_mut();
                match node.decision {
                    Decision::None => {},
                    Decision::AddForwardSearch(j) => {
                        res.added.insert(j as usize);
                        res.forward_front[0] += self.inst.p(j,0);
                        res.decision_tree = DecisionTree::add_child(&res.decision_tree, node.decision.clone());
                        for m in 1..self.inst.nb_machines() {
                            let start = max(
                                res.forward_front[m as usize - 1],
                                res.forward_front[m as usize]
                            );
                            res.forward_front[m as usize] = start + self.inst.p(j, m);
                        }
                    }
                }
            }
        }
    }
    

    pub fn add_job_forward(&self, node:&mut Node, j: JobId) -> Node {
        self.compute_lazy_part(node);  // make sure the lazy part is computed
        let node_lazypart = node.lazy_part.lazy_get();
        let lazy_part = node_lazypart.as_ref().borrow();
        let mut res = node.clone();
        // update front & idle time
        let mut front = lazy_part.forward_front[0] + self.inst.p(j,0);
        for m in 1..self.inst.nb_machines() {
            let p = self.inst.p(j, m);
            let start = max(front, lazy_part.forward_front[m as usize]);
            let current_idle = max(0, start - lazy_part.forward_front[m as usize]);
            res.idletime += current_idle;
            front = start + p;
        }
        res.decision = Decision::AddForwardSearch(j);
        res.nb_added += 1;
        let last_front;
        if node.nb_added == 0 {
            last_front = 0;
        } else {
            last_front = *lazy_part.forward_front.last().unwrap();
        }
        res.bound += front;
        res.fcost += (front-last_front)*(self.inst.nb_jobs()-node.nb_added as u32) as Time;
        res
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

