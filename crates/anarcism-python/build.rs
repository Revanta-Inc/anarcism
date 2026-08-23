fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    // CPython extension symbols resolve at load time on macOS.
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("macos") {
        println!("cargo:rustc-cdylib-link-arg=-undefined");
        println!("cargo:rustc-cdylib-link-arg=dynamic_lookup");
    }
}
