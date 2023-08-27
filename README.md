# P.H.A.J. - Probe Hybridization Analyzer in Julia
Basic filtering of potential DNA hybridization probes based on GC content, homodimer and heterodimer formation, and melting temperature.

[![Documentation](https://github.com/cookienocreams/Phaj/actions/workflows/documentation.yml/badge.svg)](https://github.com/cookienocreams/Phaj/actions/workflows/documentation.yml)

## Installation instructions

Create executable to run on local machine using the `julia` library `PackageCompiler`:

You will need to have Julia installed on your computer before starting. Julia can be installed from here: https://julialang.org/downloads/

Download the git repository using git or manually and change into the repository folder.
```bash
git clone https://github.com/cookienocreams/Phaj.git phaj
cd phaj
```
Activate the downloaded `julia` environment.
```julia
using Pkg
Pkg.activate("./")
```
The next step is to install all libraries and their dependencies.
```julia
Pkg.instantiate()
```

The last step is to create the precompiled executable. Make sure to set the correct paths for your machine.

```julia
using PackageCompiler
PackageCompiler.create_app("./", "/home/user/Phaj_app", incremental=true, precompile_execution_file="./src/Phaj.jl", include_lazy_artifacts=true)
```

There are numerous options that can be changed if desired. Use `-h` or `--help` flags to see options.

```bash
/home/user/Phaj_app/bin/Phaj --help
```

## Basic usage

The app can be run using the `Phaj` executable in the `/home/user/Phaj_app/bin` folder in a folder containing fasta to be analyzed. Note probes are expected to be in a standard fasta format with a different name for each probe.

```bash
/home/user/Phaj_app/bin/Phaj sample.fasta
```

Using the sample input fasta file, `all_probes.fa`, in `test/`, the program can be test on a dataset which contains 5000 randomly generated 80 bp potential probes. 
This can also serve as a benchmark and timing guide. The program finishes in ~ 2 minutes on a 12th Gen Intel i7-1255U, 16 GB RAM, laptop.

```bash
cd test
time /home/user/Phaj_app/bin/Phaj all_probes.fa
```
