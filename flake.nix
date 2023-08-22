{
  inputs = {
    nixpkgs.url = "nixpkgs/nixos-unstable";
    utils.url = "github:vale981/hiro-flake-utils";
    poetry2nix.url = "github:nix-community/poetry2nix";
    utils.inputs.poetry.follows = "poetry2nix";
    utils.inputs.nixpkgs.follows = "nixpkgs";
  };

  outputs = { self, utils, nixpkgs, ... }:
    (utils.lib.poetry2nixWrapper nixpkgs {
      name = "cmc_formulas";
      shellPackages = pkgs: with pkgs; [ black pyright python311Packages.jupyter pyright ];
      python = pkgs: pkgs.python311Full;
      # shellOverride = (oldAttrs: {
      #   shellHook = ''
      #               export PYTHONPATH=/home/hiro/src/hops/:$PYTHONPATH
      #               export PYTHONPATH=/home/hiro/src/hopsflow/:$PYTHONPATH
      #               export PYTHONPATH=/home/hiro/src/stocproc/:$PYTHONPATH
      #               '';
      # });
      poetryArgs = {
        projectDir = ./.;
      };
    });
}
