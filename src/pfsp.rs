use std::fs::File;
use std::io::{BufRead, BufReader, Error};
use std::vec::Vec;

pub type JobId = u32;
pub type MachineId = u32;
pub type Time = i64;

/// problem objective
pub enum ProblemObjective {
    Makespan,
    Flowtime,
}

#[derive(Debug)]
pub struct Instance {
    /// name of the instance
    name: String,
    /// number of jobs
    n: JobId,
    /// number of machines
    m: MachineId,
    /// nb_jobs / nb_machines
    job_machine_ratio: f64,
    /// processing time of job i on machine j
    p: Vec<Vec<Time>>,
    /// sum_p[m]: sum of completion times on machine m
    sum_p: Vec<Time>,
    /// sum of processing times of jobs better on the first machine
    sum_better_first: Time,
    /// sum of processing times of jobs better on the last machine
    sum_better_last: Time,
    /// simple bound using processing times
    processing_time_bound: Time,
}

impl Instance {
    /// Reads a permutation flowshop instance file
    pub fn new(filename: &str) -> Result<Instance, Error> {
        let input = File::open(filename)?;
        let buffered = BufReader::new(input);
        let mut n = 0;
        let mut m = 0;
        let mut p: Vec<Vec<Time>> = Vec::new();
        for (i, line) in buffered.lines().enumerate() {
            if i == 0 {  // first line, read n and m
                let header_line:Vec<u32> = line.unwrap()
                    .split(' ')
                    .filter(|e| !e.is_empty())
                    .map(|e| e.parse::<u32>().unwrap())
                    .collect();
                n = header_line[0];
                m = header_line[1];
            } else {  // read processing times matrix
                if ( i as MachineId ) <= m {
                    p.push(
                        line.unwrap()
                            .split(' ')
                            .filter(|e| !e.is_empty())
                            .map(|e| e.parse::<Time>().unwrap())
                            .collect()
                    );
                }
            }
        }
        let sum_p = Self::compute_sum_p(&p);
        let processing_time_bound:Time = *sum_p.iter().max().unwrap();
        let (sum_better_first,sum_better_last) = Self::compute_better_first_last(&p);
        // return the instance
        Ok(Instance {
            name: filename.to_string(),
            n,
            m,
            job_machine_ratio: n as f64 / m as f64,
            p,
            sum_p,
            sum_better_first,
            sum_better_last,
            processing_time_bound,
        })
    }

    /// creates an instance from a processing time matrix
    #[allow(dead_code)]
    pub fn from_pmj(p:Vec<Vec<Time>>) -> Self {
        let sum_p = Self::compute_sum_p(&p);
        let n = p[0].len() as JobId;
        let m = p.len() as MachineId;
        let processing_time_bound:Time = *sum_p.iter().max().unwrap();
        let (sum_better_first,sum_better_last) = Self::compute_better_first_last(&p);
        Self {
            name: "NONAME_INSTANCE".to_string(),
            n,
            m,
            job_machine_ratio: n as f64 / m as f64,
            p,
            sum_p,
            sum_better_first,
            sum_better_last,
            processing_time_bound,
        }
    }

    /// returns the number of jobs
    pub fn nb_jobs(&self) -> JobId { self.n }

    /// returns the number of machines
    pub fn nb_machines(&self) -> MachineId { self.m }

    /// returns the processing time of job j on machine m (j first, m second)
    pub fn p(&self, j:JobId, m:MachineId) -> Time { self.p[m as usize][j as usize] }

    /// returns the processing time of job j on machine m (m first, j second)
    #[allow(dead_code)]
    pub fn pmj(&self, m:MachineId, j:JobId) -> Time { self.p(j,m) }

    /// returns the sum of processing times on machine m (∑_j p(j,m))
    pub fn sum_p(&self, m:MachineId) -> Time { self.sum_p[m as usize] }

    /// job processing times lower bound
    pub fn processing_time_bound(&self) -> Time { self.processing_time_bound }

    /// true if the job j is shorter on the first machine
    pub fn better_on_first_machine(&self, j:JobId) -> bool {
        self.p(j, 0) < self.p(j, self.nb_machines()-1)
    }

    /// sum of processing times of jobs better on the first machine
    pub fn sum_better_first_machine(&self) -> Time {
        self.sum_better_first
    }

    /// sum of processing times of jobs better on the last machine
    pub fn sum_better_last_machine(&self) -> Time {
        self.sum_better_last
    }


    /* INTERNAL METHODS */
    /**
    * p[m][j]: processing time of job j on machine m
    */
    fn compute_sum_p(p: &[Vec<Time>]) -> Vec<Time> {
        let mut res = vec![0;p.len()];
        for (m,pm) in p.iter().enumerate() {
            for pmj in pm {
                res[m] += pmj;
            }
        }
        res
    }

    /// computes sums better on first/last machine
    fn compute_better_first_last(p: &[Vec<Time>]) -> (Time,Time) {
        let mut sum_first_machine:Time = 0;
        let mut sum_last_machine:Time = 0;
        let last_machine = p.len()-1;
        let n = p[0].len();
        for j in 0..n {
            if p[0][j] < p[last_machine][j] {
                sum_first_machine += p[0][j];
            } else {
                sum_last_machine += p[last_machine][j];
            }
        }
        (sum_first_machine, sum_last_machine)
    }

}

pub fn checker(inst: &Instance, objective:ProblemObjective, sol:&[JobId]) -> Time {
    let m = inst.nb_machines();
    let n = inst.nb_jobs();
    // check that every job is scheduled
    let mut is_scheduled = vec![false;n as usize];
    let mut nb_different = 0;
    for j in sol {
        if j >= &n { panic!("job {} is larger than the job number", j); }
        if !is_scheduled[*j as usize] {
            is_scheduled[*j as usize] = true;
            nb_different += 1;
        } else {
            panic!("job {} is scheduled twice", j);
        }
    }
    if nb_different != n {
        panic!("number of different jobs ({}) is different than n ({})", nb_different, n);
    }
    // compute cost
    let mut res = 0;
    let mut front = vec![0 ; m as usize];
    for j in sol {
        front[0] += inst.p(*j,0);
        for i in 1..m {
            front[i as usize] = std::cmp::max(
                front[(i-1) as usize],
                front[i as usize]
            )+inst.p(*j, i);
        }
        if let ProblemObjective::Flowtime = objective { res += front[(m-1) as usize]; };
    }
    match objective {
        ProblemObjective::Flowtime => res,
        ProblemObjective::Makespan => front[(m-1) as usize]
    }
}