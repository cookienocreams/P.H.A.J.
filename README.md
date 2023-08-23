# DNA_probe_filter
Basic filtering of potential DNA hybridization probes

## Installation instructions

Create executable to run on local machine using the `julia` library `PackageCompiler`:

You will need to have Julia installed on your computer before starting. Julia can be installed from here: https://julialang.org/downloads/

Download the git repository using git or manually and change into the repository folder.
```bash
git clone https://github.com/cookienocreams/Hybridization_probe_filter.git dna_probe_filter
cd dna_probe_filter
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
PackageCompiler.create_app("./", "/home/user/probe_filter_app", incremental=true, precompile_execution_file="./src/dna_probe_filter.jl", include_lazy_artifacts=true)
```

The app can be run using the `dna_probe_filter` executable in the `/home/user/probe_filter_app/bin` folder in a folder containing fasta to be analyzed. Note probes are expected to be in a standard fasta format with a different name for each probe.

```bash
cd fastqs
/home/user/probe_filter_app/bin/dna_probe_filter
```
