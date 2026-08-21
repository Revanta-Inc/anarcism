fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    // A CPython extension module resolves the interpreter's symbols at load
    // time, not link time. macOS refuses undefined symbols by default, so the
    // linker has to be told. `rustc-cdylib-link-arg` scopes this to this
    // crate's cdylib rather than to every target in the workspace, which keeps
    // a bare `cargo build`/`cargo clippy` working without a libpython.
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("macos") {
        println!("cargo:rustc-cdylib-link-arg=-undefined");
        println!("cargo:rustc-cdylib-link-arg=dynamic_lookup");
    }
}
