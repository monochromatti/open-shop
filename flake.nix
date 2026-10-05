{
  description = "OpenSHOP hydropower scheduling library";
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/419fe0f449b3fbe3bdd53d9840288db4509ec32e";
  outputs = { self, nixpkgs }: let
    systems = [ "aarch64-darwin" "x86_64-darwin" "aarch64-linux" "x86_64-linux" ];
  in {
    devShells = nixpkgs.lib.genAttrs systems (system: let
      pkgs = import nixpkgs { inherit system; };
    in { default = pkgs.mkShell {
      packages = [ pkgs.julia-bin (pkgs.python3.withPackages (ps: [ ps.pyyaml ])) ];
      shellHook = ''
        export JULIA_DEPOT_PATH="$PWD/.depot"
        export JULIA_NUM_THREADS=1
        export JULIA_NUM_PRECOMPILE_TASKS=2
        export OPENBLAS_NUM_THREADS=1
      '';
    }; });
  };
}
