{
  description = "Development environment for Lázer lattice library";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    utils.url = "github:numtide/flake-utils";
    # Track FLINT main directly to access newly merged features
    flint-src = {
      url = "github:flintlib/flint/main";
      flake = false;
    };
  };

  outputs = { self, nixpkgs, utils, flint-src }:
    utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };
        
        # Override the existing flint package to build from the main branch
        flint-main = pkgs.flint.overrideAttrs (oldAttrs: {
          version = "main-latest";
          src = flint-src;
          NIX_ENFORCE_NO_NATIVE = 0;

          # Ensure any patches applied to older versions are disabled 
          # if they conflict with the updated main branch code layout
          patches = []; 

          configureFlags = (oldAttrs.configureFlags or []) ++ [
            "--enable-reentrant"
          ];
        });

        sage-fast = pkgs.sage.overrideAttrs (oldAttrs: {
          doInstallCheck = false;
        });
      in
      {
        devShells.default = pkgs.mkShell {
          # Build dependencies for Lázer
          buildInputs = with pkgs; [
            flint-main # Using the modified package tracking the main branch
            gmp
            mpfr
            gnumake
            gcc
            cmake
            # sage-fast
            valgrind
          ];

          # Development tools for Vim/Neovim
          nativeBuildInputs = with pkgs; [
            clang-tools # Provides clangd (LSP)
            pkg-config
            bear        # To generate compile_commands.json
          ];

          shellHook = ''
            echo "Lázer development environment loaded with FLINT (main branch)."
            echo "Run 'bear -- make' to generate LSP metadata for Vim."
          '';
        };
      });
}
