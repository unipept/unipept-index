fn main() {
    std::process::Command::new("git")
        .args([
            "submodule",
            "update",
            "--init",
            "--recursive"
        ])
        .output()
        .expect("Failed to fetch git submodules!");
}
