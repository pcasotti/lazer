{
  description = "Development environment for Lázer lattice library";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, utils }:
    utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };
      in
      {
        devShells.default = pkgs.mkShell {
          # Build dependencies for Lázer
          buildInputs = with pkgs; [
            flint
            gmp
            mpfr
            gnumake
            gcc
            cmake
            sage
            valgrind
          ];

          # Development tools for Vim/Neovim
          nativeBuildInputs = with pkgs; [
            clang-tools # Provides clangd (LSP)
            pkg-config
            bear        # To generate compile_commands.json
          ];

          shellHook = ''
            echo "Lázer development environment loaded."
            echo "Run 'bear -- make' to generate LSP metadata for Vim."
          '';
        };
      });
}
