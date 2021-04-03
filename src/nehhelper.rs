/**
 * helper functions for the NEH algorithm
 */
use crate::pfsp::{JobId, MachineId, Instance, Time};
use std::cmp::max;


/**
 * Computes job positions given a vector of job. (i.e. compute_job_position([2,0,1]) = [1,2,0])
 */
 #[allow(dead_code)]
pub fn compute_job_positions(joblist:Vec<JobId>) -> Vec<usize> {
      let n:usize = joblist.len();
      let mut res = vec![0; n];
      for (i,e) in joblist.iter().enumerate() {
            debug_assert!((*e as usize) < n);
            res[*e as usize] = i;
      }
      return res;
}


/**
 * e, q, and, f values for the Taillard's accelerations
 */
#[derive(Debug)]
pub struct EQF {
      pub e:Vec<Vec<i64>>,
      pub q:Vec<Vec<i64>>,
      pub f:Vec<Vec<i64>>,
}

/* void calculate_e_q_f(int m,int k, int ** p_ij,int * partial_sequence,int next_job,long int ** e,long int ** q, long int ** f) */
pub fn compute_eqf(inst:&Instance, partial_sequence:&Vec<JobId>, next_job:JobId) -> EQF {
      let m:usize = inst.nb_machines() as usize;
      let k:usize = partial_sequence.len();
      let mut e = vec![vec![0;k];m as usize];
      let mut q = vec![vec![0;k+1];m as usize];
      let mut f = vec![vec![0;k+1];m as usize];
      // E COMPUTATION (earliest starting time)
      // compute earliest starting time (e)
      e[0][0] = inst.pmj(0, partial_sequence[0]);
      // compute e for the first job
      for i in 1..m {
            e[i][0] = e[(i-1) as usize][0] + inst.pmj(i as MachineId, partial_sequence[0]);
      }
      // complete for other jobs
      for j in 1..k {
            // first machine
            e[0][j] = e[0][j-1] + inst.pmj(0, partial_sequence[j]);
            // other machines
            for i in 1..m {
                  e[i][j] = max(e[i-1][j], e[i][j-1]) + inst.pmj(i as MachineId, partial_sequence[j]);
            }
      }
      // Q COMPUTATION (tail)
      // compute q for the last job
      q[m-1][k-1] = inst.pmj((m-1) as MachineId, partial_sequence[k-1]);
      for i in (0..m-1).rev() {
            q[i][k-1] = q[i+1][k-1] + inst.pmj(i as MachineId, partial_sequence[k-1]);
      }
      // compute q for the other jobs
      for j in (0..k-1).rev() {
            q[m-1][j] = q[m-1][j+1] + inst.pmj((m-1) as MachineId, partial_sequence[j]); // last machine
            // other machines
            for i in (0..m-1).rev() {
                  q[i][j] = max(q[i+1][j], q[i][j+1]) + inst.pmj(i as MachineId, partial_sequence[j])
            }
      }  
      // F COMPUTATION EARLIEST RELATIVE COMPLETION TIME (for all positions, where the job can be inserted)
      f[0][0] = inst.pmj(0,next_job);
      // first position for the rest of machines
      for i in 1..m {
            f[i][0] = f[i-1][0] + inst.pmj(i as MachineId,next_job);
      }
      // other positions
      for l in 1..k+1 {
            f[0][l] = e[0][l-1] + inst.pmj(0, next_job); // first machine
            // other machines
            for i in 1..m {
                  f[i][l] = max(e[i][l-1], f[i-1][l]) + inst.pmj(i as MachineId, next_job);
            }
      }      
      return EQF { e, q, f };
}


/**
 * computes the idle time when inserting a job
 * includes front-delays and excludes back-delays
 */
pub fn compute_idle_time(inst:&Instance, partial_sequence:&Vec<JobId>, pos:usize, eqf:&EQF) -> i64 {
      let mut res:i64 = 0;
      let m:usize = inst.nb_machines() as usize;
      let k:usize = partial_sequence.len();
      if pos < k {
            let mut fprime:Vec<i64> = vec![0;m];
            fprime[0] = eqf.f[0][pos] + inst.pmj(0, partial_sequence[pos]);
            for i in 1..m {
                  fprime[i] = max(
                        fprime[i-1], 
                        eqf.f[i][pos]+inst.pmj(i as MachineId, partial_sequence[pos])
                  );
                  res += eqf.f[i][pos]
                        - eqf.e[i][pos]
                        + inst.pmj(i as MachineId, partial_sequence[pos]) 
                        + max(0, fprime[i-1]-eqf.f[i][pos]);
            }
      } else { // inserted at the last position
            for i in 1..m {
                  res += eqf.f[i][pos] - eqf.e[i][pos-1];
            }
      }
      return res;
}


pub struct LocalState {
    pub s:Vec<JobId>,
    pub v:Time
}

/**
 * Implements an insertion neighborhood. Try to insert every task id to every other position in O(n^2)
 * This procedure uses the taillard acceleration
 * **inst**: instance
 * **s**: sequence to perturbate
 * **v**: objective value of s
*/
pub fn first_insertion_neighborhood(inst:&Instance, state:&LocalState) -> Option<LocalState> {
    for i in 0..state.s.len() { // position where to remove a job
        // create a new sub-sequence without a job
        let mut s2:Vec<JobId> = Vec::with_capacity(state.s.len());
        for j in 0..state.s.len() {
            if j != i {
                s2.push(state.s[j]);
            }
        }
        let job = state.s[i];
        // find the best possible insertion for this job
        let eqf:EQF = compute_eqf(inst, &s2, job);
        for l in 0..state.s.len() { // all possible new positions
            let mut bound:Time = eqf.f[0][l] + eqf.q[0][l];
            for i in 1..inst.nb_machines() {
                bound = max(bound, eqf.f[i as usize][l]+eqf.q[i as usize][l]);
            }
            if bound < state.v { // if a better solution found, stop there
                s2.insert(l, job);
                return Some(LocalState {
                    s: s2,
                    v: bound
                });
            }
        }
        
    }
    return None;
}


/*
 * UNIT TESTING
 */
#[cfg(test)]
mod tests {
      use super::*;

      #[test]
      fn test_compute_job_positions_a() {
            let res:Vec<usize> = compute_job_positions(vec![1,0,2]);
            assert_eq!(res[0], 1);
            assert_eq!(res[1], 0);
            assert_eq!(res[2], 2);
      }

      /**
       * 122333
       * .112..33
       * ...1112233
       */
      #[test]
      fn test_compute_eqf_e() {
            let inst:Instance = Instance::from_pmj(vec![vec![1,2,3,1],vec![2,1,2,1],vec![2,2,2,1]]);
            let partial_sequence = vec![0,1,2];  // only node 0 and 1 is added
            let res:EQF = compute_eqf(&inst, &partial_sequence, 3);
            // check that the earliest completion time of the first scheduled job is 0
            assert_eq!(res.e[0][0], 1);
            assert_eq!(res.e[1][0], 1+2);
            assert_eq!(res.e[2][0], 1+2+2);
            // second job
            assert_eq!(res.e[0][1], 3);
            assert_eq!(res.e[1][1], 4);
            assert_eq!(res.e[2][1], 7);
            // last job
            assert_eq!(res.e[0][2], 6);
            assert_eq!(res.e[1][2], 8);
            assert_eq!(res.e[2][2], 10);
      }

      #[test]
      fn test_compute_eqf_q() {
            let inst:Instance = Instance::from_pmj(vec![vec![1,2,3,1],vec![2,1,2,1],vec![2,2,2,1]]);
            let partial_sequence = vec![0,1,2];  // only node 0 and 1 is added
            let res:EQF = compute_eqf(&inst, &partial_sequence, 3);
            assert_eq!(res.q[0][0], 10);
            assert_eq!(res.q[1][0], 8);
            assert_eq!(res.q[2][0], 6);
            // second job
            assert_eq!(res.q[0][1], 9);
            assert_eq!(res.q[1][1], 5);
            assert_eq!(res.q[2][1], 4);
            // last job
            assert_eq!(res.q[0][2], 7);
            assert_eq!(res.q[1][2], 4);
            assert_eq!(res.q[2][2], 2);
      }

      #[test]
      fn test_compute_eqf_f() {
            let inst:Instance = Instance::from_pmj(vec![vec![1,2,3,1],vec![2,1,2,1],vec![2,2,2,1]]);
            let partial_sequence = vec![0,1,2];  // only node 0 and 1 is added
            let res:EQF = compute_eqf(&inst, &partial_sequence, 3);
            // first position
            assert_eq!(res.f[0][0], 1);
            assert_eq!(res.f[1][0], 2);
            assert_eq!(res.f[2][0], 3);
            //second position
            assert_eq!(res.f[0][1], 2);
            assert_eq!(res.f[1][1], 4);
            assert_eq!(res.f[2][1], 6);
            // third 
            assert_eq!(res.f[0][2], 4);
            assert_eq!(res.f[1][2], 5);
            assert_eq!(res.f[2][2], 8);
            // fourth
            assert_eq!(res.f[0][3], 7);
            assert_eq!(res.f[1][3], 9);
            assert_eq!(res.f[2][3], 11);
      }

}