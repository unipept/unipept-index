use clap::Parser;
use sa_index::{run, Arguments};

fn main() {
    let args = Arguments::parse();
    if let Err(error) = run(args) {
        eprintln!("{}", error);
        std::process::exit(1);
    };
}
