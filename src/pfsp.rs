use std::fs::File;
use std::io::{BufRead, BufReader, Error};
use std::vec::Vec;
use std::cmp::{min,max};

pub type JobId = u32;
pub type MachineId = u32;
pub type Time = i64;

#[derive(Debug)]
pub struct Instance {
    name: String,
    n: JobId,
    m: MachineId,
    /// processing time of job i on machine j
    p: Vec<Vec<Time>>,
    /// sum_p[m]: sum of completion times on machine m
    sum_p: Vec<Time>,
    /// sum of processing times of jobs better on the first machine
    sum_better_first: Time,
    /// sum of processing times of jobs better on the last machine
    sum_better_last: Time,
    /// job_spans[j] -> span of job j
    job_spans: Vec<Time>,
    /// sorted_by_span[i] -> ith job sorted by increasing span
    sorted_by_span: Vec<JobId>,
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
                    .split(" ")
                    .filter(|e| e.len() > 0)
                    .map(|e| e.parse::<u32>().unwrap())
                    .collect();
                n = header_line[0] as JobId;
                m = header_line[1] as MachineId;
            } else {  // read processing times matrix
                if ( i as MachineId ) <= m {
                    p.push(
                        line.unwrap()
                            .split(" ")
                            .filter(|e| e.len() > 0)
                            .map(|e| e.parse::<Time>().unwrap())
                            .collect()
                    );
                }
            }
        }
        let sum_p = Self::compute_sum_p(&p);
        // compute job spans
        let mut sum_first_machine:Time = 0;
        let mut sum_last_machine:Time = 0;
        let mut min_span:Time = Time::MAX;
        let mut max_span:Time = 0;
        let mut job_spans:Vec<Time> = Vec::with_capacity(n as usize);
        let last_machine = (m-1) as usize;
        for j in 0..n {
            let mut tmp:Time = 0;
            if p[0][j as usize] < p[last_machine][j as usize] {
                sum_first_machine += p[0][j as usize];
            } else {
                sum_last_machine += p[last_machine][j as usize];
            }
            for i in 0..m {
                tmp += p[i as usize][j as usize];
            }
            job_spans.push(tmp);
            min_span = min(min_span, tmp);
            max_span = max(max_span, tmp);
        }
        println!("min/max span size: {}\t{}", min_span, max_span);
        let mut sorted_jobs_by_span:Vec<JobId> = Vec::with_capacity(n as usize);
        for j in 0..n {
            sorted_jobs_by_span.push(j);
        }
        sorted_jobs_by_span.sort_by_key(|j| job_spans[*j as usize]);
        sorted_jobs_by_span.reverse();
        println!("{}\t{}\t{}",
            job_spans[sorted_jobs_by_span[0] as usize],
            job_spans[sorted_jobs_by_span[1] as usize],
            job_spans[sorted_jobs_by_span[2] as usize]
        );
        // return the instance
        Ok(Instance {
            name: filename.to_string(),
            n,
            m,
            p,
            sum_p,
            sum_better_first: sum_first_machine,
            sum_better_last: sum_last_machine,
            job_spans,
            sorted_by_span: sorted_jobs_by_span,
        })
    }

    /// creates an instance from a processing time matrix
    #[allow(dead_code)]
    pub fn from_pmj(p:Vec<Vec<Time>>) -> Self {
        let sum_p = Self::compute_sum_p(&p);
        Self {
            name: "NONAME_INSTANCE".to_string(),
            n: p[0].len() as JobId,
            m: p.len() as MachineId,
            p,
            sum_p,
            sum_better_first: 0, // TODO
            sum_better_last: 0, // TODO
            job_spans: Vec::new(), // TODO
            sorted_by_span: Vec::new(), // TODO
        }
    }

    pub fn nb_jobs(&self) -> JobId { self.n }

    pub fn nb_machines(&self) -> MachineId { self.m }

    pub fn p(&self, j:JobId, m:MachineId) -> Time { self.p[m as usize][j as usize] }

    #[allow(dead_code)]
    pub fn pmj(&self, m:MachineId, j:JobId) -> Time { self.p(j,m) }

    pub fn sum_p(&self, m:MachineId) -> Time { self.sum_p[m as usize] }

    pub fn better_on_first_machine(&self, j:JobId) -> bool {
        self.p(j, 0) < self.p(j, self.nb_machines()-1)
    }

    pub fn sum_better_first_machine(&self) -> Time {
        self.sum_better_first
    }

    pub fn sum_better_last_machine(&self) -> Time {
        self.sum_better_last
    }

    // pub fn ith_sorted_by_span(&self, i:usize) -> JobId {
    //     self.sorted_by_span[i]
    // }

    // pub fn job_span(&self, j:JobId) -> Time {
    //     self.job_spans[j as usize]
    // }

    /* INTERNAL METHODS */
    /**
    * p[m][j]: processing time of job j on machine m
    */
    fn compute_sum_p(p: &Vec<Vec<Time>>) -> Vec<Time> {
        let mut res = vec![0;p.len()];
        for m in 0..p.len() {
            for j in 0..p[m].len() {
                res[m] += p[m][j];
            }
        }
        res
    }

}

