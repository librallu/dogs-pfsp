use std::fs::File;
use std::io::Write;
use std::rc::Rc;
use std::time::SystemTime;
use structopt::StructOpt;

use dogs::searchspace::SearchSpace;
use dogs::treesearch::algo::iterative_beamsearch::IterativeBeamSearch;
use dogs::treesearch::decorators::stats::StatTsDecorator;
// use dogs::treesearch::decorators::pruning::PruningDecorator;
// use dogs::treesearch::decorators::bounding::BoundingDecorator;

// test logger
use dogs::metriclogger::MetricLogger;

mod pfsp;
// mod fbwamakespan;
// mod fbwamakespan2;
// mod ffg3flowtime;

mod fbsearch;

// mod nehhelper;
// mod nehmakespan;
// use nehmakespan::{NEHSearch,TiebreakingScheme};

use fbsearch::{BidirectionalScheme, BranchingScheme, FBSearch, Guide, Objective};

#[derive(Debug, StructOpt)]
struct Cli {
    /// path to instance file
    instance: String,
    /// time for the algorithm to run
    t: f32,
    /// objective to minimize (values: [makespan, flowtime])
    obj: String,
    /// guide function (values: [bound, idle, alpha, walpha, wfrontalpha, gap])
    guide: String,
    /// branching scheme (values: [forward, bi_min, bi_alt])
    branchingscheme: String,
    /// output prefix (will write a few files starting by this prefix)
    output: String,
}

fn s2obj(s: &String) -> Option<Objective> {
    if s == "makespan" {
        return Some(Objective::Makespan);
    } else if s == "flowtime" {
        return Some(Objective::Flowtime);
    } else {
        return None;
    }
}

fn s2guide(s: &String) -> Option<Guide> {
    if s == "bound" {
        return Some(Guide::Bound);
    } else if s == "idle" {
        return Some(Guide::IdleTime);
    } else if s == "alpha" {
        return Some(Guide::BoundAlphaIdle);
    } else if s == "walpha" {
        return Some(Guide::BoundW2AlphaIdle);
    } else if s == "wfrontalpha" {
        return Some(Guide::BoundWAlphaIdle);
    } else if s == "gap" {
        return Some(Guide::BoundGapIdle);
    } else {
        return None;
    }
}

fn s2branching(s: &String) -> Option<BranchingScheme> {
    if s == "forward" {
        return Some(BranchingScheme::Forward);
    } else if s == "bi_min" {
        return Some(BranchingScheme::Bidirectional(
            BidirectionalScheme::MinChildren,
        ));
    } else if s == "bi_alt" {
        return Some(BranchingScheme::Bidirectional(
            BidirectionalScheme::Alternate,
        ));
    } else {
        return None;
    }
}

fn main() {
    // start computing time
    let start_time = SystemTime::now();

    // read command line arguments
    let args = Cli::from_args();

    // create logger
    let logger = Rc::new(MetricLogger::new());

    // create search space
    let mut space = StatTsDecorator::new(FBSearch::new(
        &args.instance,
        s2obj(&args.obj).unwrap(),
        s2guide(&args.guide).unwrap(),
        s2branching(&args.branchingscheme).unwrap(),
        Some((&args.output).clone() + ".sol.txt"),
    ))
    .bind_logger(Rc::downgrade(&logger));

    // create search space with bounding decorator
    // let mut space = BoundingDecorator::new(
    //     StatTsDecorator::new(FBSearch::new(
    //         &args.instance,
    //         s2obj(&args.obj).unwrap(),
    //         s2guide(&args.guide).unwrap(),
    //         s2branching(&args.branchingscheme).unwrap(),
    //         Some((&args.output).clone() + ".sol.txt"),
    //     ))
    //     .bind_logger(Rc::downgrade(&logger))
    // ).bind_logger(Rc::downgrade(&logger));

    // let mut space = StatTsDecorator::new(
    //     NEHSearch::new(
    //         &args.instance,
    //         TiebreakingScheme::FF
    //     )
    // ).bind_logger(Rc::downgrade(&logger));

    // create the search algorithm
    let mut ts = IterativeBeamSearch::new(&mut space, 1, 2.).bind_logger(Rc::downgrade(&logger));

    // run the search algorithm and display the results afterwards
    logger.display_headers();
    ts.run(|_| start_time.elapsed().unwrap().as_secs_f32() < args.t);
    // if optimal solution found
    if !ts.is_heuristic_pruning_done() {
        let filename: String = (&args.output).clone() + ".optimal.txt";
        let mut file = match File::create(&filename) {
            Err(why) => panic!("couldn't create {}: {}", filename, why),
            Ok(file) => file,
        };
        writeln!(&mut file, "Y").unwrap();
    }
    space.display_statistics();

    // export statistics to file
    MetricLogger::export_json_stats(&space, (&args.output).clone() + ".stats.json");
    let perfprofile_filename = (&args.output).clone() + ".perfprofile.json";
    let perfprofile_name: String = format!("{} {}", &args.guide, &args.branchingscheme,);
    space.export_performance_profile(perfprofile_filename.as_str(), perfprofile_name.as_str());
}
