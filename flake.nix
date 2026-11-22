{
  description = "OpenSHOP hydropower scheduling library";
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/419fe0f449b3fbe3bdd53d9840288db4509ec32e";
  inputs.rust-overlay = {
    url = "github:oxalica/rust-overlay";
    inputs.nixpkgs.follows = "nixpkgs";
  };
  inputs.crane.url = "github:ipetkov/crane";
  outputs = { self, nixpkgs, rust-overlay, crane }: let
    systems = [ "aarch64-darwin" "x86_64-darwin" "aarch64-linux" "x86_64-linux" ];
    rust = system: let
      pkgs = import nixpkgs { inherit system; overlays = [ rust-overlay.overlays.default ]; };
      toolchain = pkgs.rust-bin.fromRustupToolchainFile ./rust-toolchain.toml;
      craneLib = (crane.mkLib pkgs).overrideToolchain toolchain;
      native = {
        nativeBuildInputs = [ pkgs.pkg-config pkgs.rustPlatform.bindgenHook ];
        buildInputs = [ pkgs.ipopt pkgs.scipopt-scip ];
        SCIPOPTDIR = "${pkgs.scipopt-scip}";
      };
    in { inherit pkgs toolchain craneLib native; };
  in {
    devShells = nixpkgs.lib.genAttrs systems (system: let
      pkgs = import nixpkgs { inherit system; };
    in {
    paper = pkgs.mkShell { packages = [ pkgs.typst pkgs.gnumake ]; };
    rust = let r = rust system; in r.pkgs.mkShell (r.native // {
      packages = [ r.toolchain r.pkgs.cargo-nextest r.pkgs.maturin
        (r.pkgs.python3.withPackages (ps: [ ps.numpy ps.pytest ps.pip ])) ]
        ++ r.pkgs.lib.optionals r.pkgs.stdenv.hostPlatform.isLinux [ r.pkgs.patchelf ];
      RAYON_NUM_THREADS = "1";
      OMP_NUM_THREADS = "1";
      OPENBLAS_NUM_THREADS = "1";
      RUSTUP_TOOLCHAIN = "${r.toolchain}";
    } // r.pkgs.lib.optionalAttrs r.pkgs.stdenv.hostPlatform.isDarwin {
      # Maturin follows transitive solver dependencies when repairing wheels.
      DYLD_LIBRARY_PATH = r.pkgs.lib.makeLibraryPath [ r.pkgs.gfortran.cc.lib ];
    });
    default = pkgs.mkShell {
      packages = [ pkgs.julia-bin (pkgs.python3.withPackages (ps: [ ps.pyyaml ])) ];
      shellHook = ''
        export JULIA_DEPOT_PATH="$PWD/.depot"
        export JULIA_NUM_THREADS=1
        export JULIA_NUM_PRECOMPILE_TASKS=2
        export OPENBLAS_NUM_THREADS=1
      '';
    }; });
    packages = nixpkgs.lib.genAttrs systems (system: let
      r = rust system;
      src = r.craneLib.cleanCargoSource ./.;
      common = r.native // {
        inherit src;
        pname = "openshop";
        version = "0.9.0";
        strictDeps = true;
        cargoExtraArgs = "-p openshop";
      };
      cargoArtifacts = r.craneLib.buildDepsOnly common;
    in {
      openshop = r.craneLib.buildPackage (common // { inherit cargoArtifacts; });
    });
  };
}
