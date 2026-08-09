use std::path::Path;
use std::process::Command;

fn main() {
    let hash = Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .ok()
        .filter(|o| o.status.success())
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .unwrap_or_else(|| "unknown".to_string());

    println!("cargo:rustc-env=GIT_COMMIT_HASH={}", hash);

    // `rerun-if-changed` paths are resolved relative to the *package* directory, not the
    // workspace root — so a bare ".git/HEAD" points at "sa-benchmarks/.git/HEAD", which does
    // not exist, and cargo then re-runs this script (and rebuilds the crate) on every
    // invocation. The repo's git dir is one level up.
    //
    // Both paths are existence-checked: emitting a watch on a missing path is what causes the
    // permanent rebuild, and the crate must still build from a source tarball with no .git.
    for path in ["../.git/HEAD", "../.git/refs/heads"] {
        if Path::new(path).exists() {
            println!("cargo:rerun-if-changed={}", path);
        }
    }
}
