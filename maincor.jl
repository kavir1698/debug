using Distributed
using Pkg
Pkg.activate(".")
addprocs(2)
@everywhere using EvoDynamics
@everywhere using Agents
@everywhere using Statistics
# using EvoDynamics
# using Agents
# using Statistics
using JLD2
using FileIO
using Distributions
@everywhere include("data_collection_functions.jl")
# include("data_collection_functions.jl")

pf = "parameters_frame.jl"
paramdir = "parameter_files_correlated/"
results_dir = "sim_outputs_correlated/"

migration_thresholds = [0.0, 0.1, 1.0, 2.0, 5.0, 10.0]
biotic_variances = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0]


pmat = trues(3, 21)  # the original pleiotropy_matrix that we will randomize a little bit
nchanges_mean = 4
pmat_variations = 10
nreplicates = 3


"""
A variation of the create single_parameter_file function that adds random 1's to the pleiotropy_matrix

## Parameters

* nchanges_mean: mean of a Poisson distribution for the number of 0s in the pleiotropy_matrix that will change to 1.
"""
function create_single_parameter_file_add_noise_to_pleiotropy2(pf, pmat, pmatnum, migration_rate, biotic_coeff, paramdir=paramdir, nchanges_mean=2)
  if !isdir(paramdir)
    mkdir(paramdir)
  end
  outfile = joinpath(paramdir, "params_migration_rate=$(migration_rate)_biotic_coeff=$(biotic_coeff)_fixed_interaction_mat_$(pmatnum).jl")
  input = readlines(pf)
  # Update migration threshold
  lines = findall(x -> startswith(strip(x), ":migration_threshold"), input)
  input[lines] .= ":migration_threshold => $(migration_rate),"

  # Update biotic variance
  lines = findall(x -> startswith(strip(x), ":biotic_variance"), input)
  input[lines] .= ":biotic_variance => $(biotic_coeff),"

  # Add some ones to the pleiotropy_matrix
  lines = findall(x -> startswith(strip(x), ":pleiotropy_matrix"), input)
  nchanges_dist = Poisson(nchanges_mean)
  nchanges = rand(nchanges_dist)
  pmat2 = deepcopy(pmat)
  pmat2[rand(findall(x -> x == true, pmat), nchanges)] .= false
  input[lines] .= ":pleiotropy_matrix => $(pmat2),"

  open(outfile, "w") do ff
    for line in input
      println(ff, line)
    end
  end
  return outfile
end

function create_parameter_combinations2(pf, pmat, migration_rates, biotic_coeffs, paramdir=paramdir, nchanges_mean=2)
  for mrate in migration_rates
    for biocoef in biotic_coeffs
      for pmatnum in 1:pmat_variations
        create_single_parameter_file_add_noise_to_pleiotropy2(pf, pmat, pmatnum, mrate, biocoef, paramdir, nchanges_mean)
      end
    end
  end
end

create_parameter_combinations2(pf, pmat, migration_thresholds, biotic_variances, paramdir, nchanges_mean)

all_parameter_files = readdir(paramdir)

for f in all_parameter_files
  if !isdir(results_dir)
    mkdir(results_dir)
  end
  param_file = joinpath(paramdir, f)
  adata, mdata, models = runmodel(param_file, adata=nothing, mdata=[EvoDynamics.mean_fitness_per_species, EvoDynamics.species_N, mean_espistasis_matrix_per_species], parallel=false, when_model=0:10:500, showprogress=true, offline_run=false, writing_interval=10, mdata_filename=joinpath(results_dir, "$(f[1:end-3]).csv"))
  save(joinpath(results_dir, "$(f[1:end-3]).jld2"), "results", mdata)
  break
end
