//! useful additional warnings if docs are missing, or crates imported but unused, etc.
#![warn(missing_debug_implementations)]
#![warn(missing_docs)]
#![warn(trivial_casts, trivial_numeric_casts)]
#![warn(unsafe_code)]
#![warn(unused_extern_crates)]
#![warn(variant_size_differences)]

//! not sure if already by default in clippy
#![warn(clippy::similar_names)]
#![warn(clippy::shadow_unrelated)]
#![warn(clippy::shadow_same)]
#![warn(clippy::shadow_reuse)]

//! checks integer arithmetic in the project
// #![warn(clippy::integer_arithmetic)]

//! these flags can be useful, but will indicate normal behavior
// #![warn(clippy::use_debug)]
// #![warn(clippy::print_stdout)]
// #![warn(clippy::cast_possible_truncation)]
// #![warn(clippy::cast_possible_wrap)]
// #![warn(clippy::cast_precision_loss)]
// #![warn(clippy::cast_sign_loss)]

use std::cell::RefCell;
use std::fs::File;
use std::io::Write;
use std::rc::Rc;

use ordered_float::OrderedFloat;

#[macro_use]
extern crate clap;
use clap::App;

use serde_json::json;

use dogs::metric_logger::MetricLogger;
use dogs::search_algorithm::{SearchAlgorithm, TimeStoppingCriterion};
use dogs::search_space::{SearchSpace,GuidedSpace,TotalNeighborGeneration};
use dogs::tree_search::algo::beam_search::create_iterative_beam_search;
use dogs::tree_search::decorators::pruning::PruningDecorator;
use dogs::tree_search::decorators::stats::StatTsDecorator;
// use dogs::treesearch::decorators::bounding::BoundingDecorator;


mod pfsp;
mod fflowtime;
mod fbmakespan;
mod nehhelper;
mod nehmakespan;
// mod expandmakespan;

fn run_search_tree<T,N>(space:T, t:f32, perf_file: Option<String>, inst_name: &str, algo_name:String) where
T:SearchSpace<N,pfsp::Time>+GuidedSpace<N,OrderedFloat<f64>>+TotalNeighborGeneration<N>,
N:Clone {
    // create logger and stopping criterion
    let logger = Rc::new(MetricLogger::default());
    let stopping_criterion = TimeStoppingCriterion::new(t);
    // create search space
    let space = Rc::new(RefCell::new(
        StatTsDecorator::new(
            // BoundingDecorator::new(
                space
            // ).bind_logger(Rc::downgrade(&logger))
        )
        .bind_logger(Rc::downgrade(&logger)),
    ));
    // create the search algorithm
    let mut ts =
        create_iterative_beam_search(space.clone(), 1., 2.)
        .bind_logger(Rc::downgrade(&logger));
    logger.display_headers();
    ts.run(stopping_criterion);
    // display the results afterwards
    space.borrow_mut().display_statistics();
    // export statistics to file
    if let Some(f) = perf_file {
        let mut res = json!({});
        res["inst"]["name"] = json!(inst_name);
        res["IBS"]["algo_name"] = json!(algo_name);
        space.borrow().json_statistics(&mut res["search_space"]);
        ts.json_statistics(&mut res["IBS"]);
        let mut file = match File::create(f.as_str()) {
            Err(why) => panic!("couldn't create {}: {}", f, why),
            Ok(file) => file
        };
        if let Err(why) = file.write(serde_json::to_string(&res).unwrap().as_bytes()) {
            panic!("couldn't write: {}",why)
        };
    }
}

fn main() {
    let yaml = load_yaml!("main_args.yml");
    let main_args = App::from_yaml(yaml).get_matches();
    let inst_filename = main_args.value_of("instance").unwrap();
    let t:f32 = main_args.value_of("time").unwrap().parse::<f32>().unwrap();
    let sol_file: Option<String> = match main_args.value_of("solution") {
        None => None,
        Some(e) => {
            println!("======================================");
            println!("\t printing solutions in: {}", e);
            Some(e.to_string())
        }
    };
    let perf_file: Option<String> = match main_args.value_of("perf") {
        None => None,
        Some(e) => {
            println!("\t printing perfs in: {}\n", e);
            Some(e.to_string())
        }
    };
    // construct the search algorithm

    // construct the search space
    if let Some(branching_args) = main_args.subcommand_matches("f_flowtime") {
        let guide = match branching_args.value_of("guide").unwrap() {
            "bound" => fflowtime::Guide::Bound,
            "idle" => fflowtime::Guide::Idle,
            "alpha" => fflowtime::Guide::Alpha,
            _ => {
                panic!("unknown guide (possible values: [bound, idle, alpha]")
            }
        };
        println!(" ============== Forward Flowtime(g={:?}) ==============\n", guide);
        run_search_tree(
            PruningDecorator::new(
                fflowtime::ForwardSearch::new(
                    inst_filename,
                    guide.clone(),
                    sol_file.clone(),
                )
            )
            , t, perf_file.clone(), inst_filename, format!("IBS_forward_{:?}", guide)
        );
    }
    if let Some(branching_args) = main_args.subcommand_matches("neh_makespan") {
        let guide = match branching_args.value_of("guide").unwrap() {
            "bound" => nehmakespan::Guide::Bound,
            "ff"    => nehmakespan::Guide::FF,
            "alpha" => nehmakespan::Guide::Alpha,
            _ => {
                panic!("unknown guide (possible values: [bound, ff, alpha]")
            }
        };
        println!(" ============== NEH makespan(g={:?}) ==============\n", guide);
        run_search_tree(
            PruningDecorator::new(nehmakespan::NEHSearch::new(
                inst_filename,
                guide.clone(),
                // sol_file.clone(),
            )), t, perf_file.clone(), inst_filename, format!("IBS_NEH_{:?}", guide)
        );
    }
    if let Some(branching_args) = main_args.subcommand_matches("fb_makespan") {
        let guide = match branching_args.value_of("guide").unwrap() {
            "bound" =>  fbmakespan::Guide::Bound,
            "idle" =>   fbmakespan::Guide::Idle,
            "alpha" =>  fbmakespan::Guide::Alpha,
            "walpha" => fbmakespan::Guide::Walpha,
            "gap" => fbmakespan::Guide::Gap,
            "combined" => fbmakespan::Guide::Combined,
            _ => {
                panic!("unknown guide (possible values: [bound, idle, alpha, walpha, gap]")
            }
        };
        let branching = match branching_args.value_of("branching").unwrap() {
            "forward" =>  fbmakespan::BranchingScheme::Forward,
            "bi_alt" =>   fbmakespan::BranchingScheme::Bidirectional(fbmakespan::BidirectionalScheme::Alternate),
            "bi_min" =>   fbmakespan::BranchingScheme::Bidirectional(fbmakespan::BidirectionalScheme::MinChildren),
            _ => {
                panic!("unknown branching (possible values: [forward, bi_alt, bi_min]")
            }
        };
        let algo_name:String = format!("IBS_{:?}_{:?}", branching, guide);
        run_search_tree(
            fbmakespan::FBMakespan::new(
                inst_filename,
                guide,
                branching,
                false,
                sol_file,
            ), t, perf_file, inst_filename, algo_name
        );
    }

}
