use clap::{ArgGroup, Args};
use num_cpus;
use std::env;

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

    /// Cap cores using the min value from env vars NSLOTS or IFX_LOCAL_PROCS
    #[clap(short = 'C', long)]
    pub cap_cores_using_env: bool,
}

pub fn num_procs_process(args: NumProcsArgs) -> Result<(), std::io::Error> {
    let mut cores = if args.physical {
        num_cpus::get_physical()
    } else {
        num_cpus::get()
    };

    if args.cap_cores_using_env
        && let Some(n) = ["NSLOTS", "IFX_LOCAL_PROCS"]
            .into_iter()
            .find_map(|v| env::var(v).ok().and_then(|v| v.parse::<usize>().ok()))
    {
        if n > 0 {
            cores = cores.min(n);
        } else {
            eprintln!("IRMA-core WARNING! The requested core cap {n} is not valid, ignoring and using {cores} cores.");
        }
    }

    if let Ok(s) = env::var("LOCAL_PROCS_OVERRIDE")
        && let Ok(n) = s.parse::<usize>()
    {
        if n > 0 {
            cores = n;
        } else {
            eprintln!(
                "IRMA-core WARNING! The requested 'LOCAL_PROCS_OVERRIDE={s}' is not valid, ignoring and using {cores} cores."
            );
        }
    }

    if args.include_half {
        let half_cores = (cores / 2).max(1);
        println!("{cores} {half_cores}");
    } else {
        println!("{cores}");
    }

    Ok(())
}
