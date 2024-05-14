use clap::Parser;
use sa_index::{Arguments, run};

fn main() {
    let args = Arguments::parse();
    if let Err(error) = run(args) {
        eprintln!("{}", error);
        std::process::exit(1);
    };
}
