use clap::{ArgGroup, Args};
use num_cpus;

/// Arguments for specifying the number of processor cores to use
#[derive(Args, Debug)]
#[clap(group = ArgGroup::new("cores_type").required(true))]
pub struct NumProcsArgs {
    /// Use physical core count
    #[clap(short = 'P', long, group = "cores_type")]
    pub physical: bool,

    /// Use logical core count (includes hyperthreaded cores)
    #[clap(short = 'L', long, group = "cores_type")]
    pub logical: bool,

    /// Also include half the number of cores in the output
    #[clap(short = 'H', long)]
    pub include_half: bool,
}

pub fn num_procs_process(args: NumProcsArgs) -> Result<(), std::io::Error> {
    let cores = if args.physical {
        num_cpus::get_physical()
    } else {
        num_cpus::get()
    };

    if args.include_half {
        let half_cores = cores / 2;
        println!("{cores} {half_cores}");
    } else {
        println!("{cores}");
    }

    Ok(())
}
