use std::fs::File;
use std::io::{BufRead, BufReader, Error};
use std::vec::Vec;

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
        // return the instance
        Ok(Instance {
            name: filename.to_string(),
            n,
            m,
            p,
            sum_p,
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
        }
    }

    pub fn nb_jobs(&self) -> JobId { self.n }

    pub fn nb_machines(&self) -> MachineId { self.m }

    pub fn p(&self, i:JobId, j:MachineId) -> Time { self.p[j as usize][i as usize] }

    #[allow(dead_code)]
    pub fn pmj(&self, m:MachineId, j:JobId) -> Time { self.p(j,m) }

    pub fn sum_p(&self, m:MachineId) -> Time { self.sum_p[m as usize] }


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

