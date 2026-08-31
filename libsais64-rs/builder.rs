use std::{
    env,
    error::Error,
    fmt::{Display, Formatter},
    path::{Path, PathBuf},
    process::{Command, ExitStatus}
};

/// Custom error for compilation of the C library
#[derive(Debug)]
struct CompileError<'a> {
    command: &'a str,
    exit_code: Option<i32>
}

impl Display for CompileError<'_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let end_text = if let Some(code) = self.exit_code {
            format!("with exit code {}", code)
        } else {
            "without exit code".to_string()
        };
        let text = format!("Command with name `{}` failed {}", self.command, end_text);
        write!(f, "{}", text)
    }
}

impl Error for CompileError<'_> {}

/// Handles the exit statuses of the executed bash commands
///
/// # Arguments
/// * `name` - Name of the executed bash command
/// * `exit_states` - The exit status of the executed bash command
///
/// # Returns
///
/// Returns () if the exit status was success
///
/// # Errors
///
/// Returns a CompileError if the command failed
fn exit_status_to_result(name: &str, exit_status: ExitStatus) -> Result<(), CompileError<'_>> {
    match exit_status.success() {
        true => Ok(()),
        false => Err(CompileError { command: name, exit_code: exit_status.code() })
    }
}

/// Where the C library comes from, and exactly which commit of it.
///
/// Pinned to a hash rather than tracking the default branch, so that a given Rust commit always
/// builds the same binary and no change in the C repository can reach a build — or break CI —
/// without a commit here recording it.
///
/// Bumping is deliberate: change the hash below and commit it, so the bump is reviewable and
/// bisectable like any other dependency change.
const LIBSAIS_REPOSITORY: &str = "https://github.com/unipept/libsais-packed.git";
const LIBSAIS_COMMIT: &str = "a2b99260a4c41ad101d3388fcca1d4c78fc4a0c3";
const LIBSAIS_DIRECTORY: &str = "libsais-packed";

/// Whether the working copy is already the pinned commit.
///
/// Any failure — no directory, not a repository, a git that will not run — is answered with
/// `false`, which just means the checkout below happens.
fn is_at_pinned_commit() -> bool {
    Command::new("git")
        .args(["-C", LIBSAIS_DIRECTORY, "rev-parse", "HEAD"])
        .output()
        .map(|out| out.status.success() && String::from_utf8_lossy(&out.stdout).trim() == LIBSAIS_COMMIT)
        .unwrap_or(false)
}

/// Whether the checkout has uncommitted changes to tracked files.
///
/// Untracked files are ignored: cmake and make leave their output in this directory, so it is never
/// clean by that measure.
fn has_local_changes() -> bool {
    Command::new("git")
        .args(["-C", LIBSAIS_DIRECTORY, "status", "--porcelain", "--untracked-files=no"])
        .output()
        .map(|out| out.status.success() && !out.stdout.is_empty())
        .unwrap_or(false)
}

/// Puts [`LIBSAIS_COMMIT`] of the C library in `libsais-packed/`.
///
/// A shallow fetch of the single commit by hash, which costs what a `--depth=1` clone costs; GitHub
/// serves any commit reachable from a branch this way. A directory already at that commit is left
/// alone, so an incremental build skips the fetch entirely — `main` still reconfigures with cmake
/// and runs make, but both then find their work already done.
///
/// # Errors
///
/// Any of the git invocations failing.
fn fetch_libsais() -> Result<(), Box<dyn Error>> {
    if is_at_pinned_commit() {
        // Edits to the C source are the one reason to be here by hand, and the re-checkout below
        // starts with `rm -rf`, so they are reported rather than deleted. The build is no longer
        // the pinned one, and saying so is the most this can do without eating someone's work.
        if has_local_changes() {
            println!(
                "cargo:warning=libsais64-rs: {}/ is at the pinned commit but has uncommitted changes; \
                 building those instead of the pinned source. Delete the directory for a clean checkout.",
                LIBSAIS_DIRECTORY
            );
        }
        return Ok(());
    }

    // Whatever is there is not what is wanted; if removing fails, it is because the folder did not
    // exist, and we can ignore that.
    Command::new("rm").args(["-rf", LIBSAIS_DIRECTORY]).status().unwrap_or_default();

    exit_status_to_result("git init", Command::new("git").args(["init", "-q", LIBSAIS_DIRECTORY]).status()?)?;
    exit_status_to_result(
        "git remote add",
        Command::new("git")
            .args(["-C", LIBSAIS_DIRECTORY, "remote", "add", "origin", LIBSAIS_REPOSITORY])
            .status()?
    )?;
    exit_status_to_result(
        "git fetch",
        Command::new("git")
            .args(["-C", LIBSAIS_DIRECTORY, "fetch", "-q", "--depth=1", "origin", LIBSAIS_COMMIT])
            .status()?
    )?;
    exit_status_to_result(
        "git checkout",
        Command::new("git").args(["-C", LIBSAIS_DIRECTORY, "checkout", "-q", "FETCH_HEAD"]).status()?
    )?;

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    // fetch the pinned revision of the c library
    fetch_libsais()?;

    // compile the c library
    Command::new("rm").args(["-f", "libsais-packed/CMakeCache.txt"]).status().unwrap_or_default(); // if removing fails, it is since the cmake cache did not exist, we just can ignore it
    exit_status_to_result(
        "cmake",
        Command::new("cmake")
            .args(["-DCMAKE_BUILD_TYPE=\"Release\"", "libsais-packed", "-Blibsais-packed"])
            .status()?
    )?;
    exit_status_to_result("make", Command::new("make").args(["-C", "libsais-packed"]).status()?)?;

    // link the c libsais-packed library to rust
    let dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    println!("cargo:rustc-link-search=native={}", Path::new(&dir).join("libsais-packed").display());
    println!("cargo:rustc-link-lib=static=libsais");

    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header("libsais-wrapper.h")
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        // Finish the builder and generate the bindings.
        .generate()?;

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    let out_path = PathBuf::from(env::var("OUT_DIR")?);
    bindings.write_to_file(out_path.join("bindings.rs"))?;

    Ok(())
}
