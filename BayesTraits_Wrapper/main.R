# Title: ???
# Author: [Name]
# Date Created: ???
# Description: ???

# Load libraries
library(tidyverse)
library(phytools)


# GET CURRENT FOLDER & SET PATHS ----------------------------------------------------
path = box::file()
path_data_source = box::file("data_source")
path_data_working = box::file("data_working")
path_results = box::file("results")
path_scripts = box::file("modules") 


# SET JOB PARAMETERS ----------------------------------------------------
job_name = "Test_1_VR"

input_data = file.path(path_data_source, "ToyTreeData.tsv")
input_tree = file.path(path_data_source, "ToyTree.tre")

working_data = file.path(path_data_working, "data_trait_1.tsv")
working_tree = file.path(path_data_working, "tree.tree")

# Don't edit among runs
command_file = paste0(job_name, ".BTsettings.txt")
log_file = paste0(job_name, ".Log.txt")
stones_file = paste0(job_name, ".Stones.txt")


# LOAD & PROCESS DATA ----------------------------------------------------
# Load data
df = read.delim(input_data, header=FALSE, sep = "\t", row.names = 1)

# Process data
df_1 = df[1]

# Save processed data
write.table(df_1, file=working_data, sep = "\t", row.names=TRUE, col.names=FALSE, quote=FALSE)


# LOAD & PROCESS TREE ----------------------------------------------------
# Note: The code a mulptiphylo object (set of trees)
tree = read.nexus(input_tree)

# Remove node labels and get the tree in nexus format
tree = lapply(tree, function(t) { t$node.label = NULL; t })
class(tree) = "multiPhylo"

# Save the modified tree
write.nexus(tree, file=working_tree, translate=TRUE)


# BAYESTRAITS JOB SETTINGS ----------------------------------------------------
### BayesTraits Main Options
# 1 MultiState Multistate
# 2 Discrete: Independent Two binary traits
# 3 Discrete: Dependent Two binary traits
# 4 Continuous: Random Walk Continuous traits
# 5 Continuous: Directional Continuous traits
# 6 Continuous: Regression Two or more continuous traits
# 7 Independent Contrasts Continuous traits
# 8 Independent Contrasts: Correlation Two or more continuous traits
# 9 Independent Contrasts: Regression Two or more continuous traits
# 10 Discrete: Covarion Two binary traits
# 11 Discrete: Heterogeneous Two binary traits
# 12 Geographical Two traits longitude and latitude

# Each element is one line you'd normally type at the BayesTraits prompt
# NOTE: BayesTraits' own prompt parser can't handle a path with spaces (e.g. our
# OneDrive path), so "logfile" must get a bare name — never a full path. We instead
# setwd() into path_results before running so the bare name lands in the right place.
commands = c(
  "7",                           # Model
  "2",                           # 1:ML  2:MCMC
  "prior Alpha-1 normal 10 5",
  "prior Sigma-1 normal 10 5",
  "vr",
  "Iterations 1000000",
  "Burnin 1000",
  "Stones 100 1000",
  paste0("logfile ", job_name), # bare name only — see NOTE above
  "Run"                # must be last — starts the analysis
)

# Run BayesTraits
# First ensure that output is stored in the results folder, then go back
setwd(path_results)
out = system2(
  "BayesTraitsV5",              # found via PATH; use exact exe name if different
  args   = c(shQuote(working_tree), shQuote(working_data)), # shQuote handles spaces (e.g. OneDrive path)
  input  = commands,            # piped to stdin, replacing interactive typing
  stdout = TRUE,
  stderr = TRUE
)
setwd(path)


# PROCESS BayesTraits RESULTS ----------------------------------------------------
# Loads read_bt_log(), plot_single_log(), plot_multi_run_summary(), etc.
source(file.path(path_scripts, "bt_logPlotting.R"))

# Process Results:
# 1. Density Plots
# 2. Trace Plots
# 3. Summary table for parameters
result_log_file = file.path(path_results, paste0(job_name, ".Log.txt"))
plot_single_log(result_log_file, path_results)


# PROCESS VARIABLE RESULTS ----------------------------------------------------
# If taxa have been dropped in VR analysis, they must be be dropped from the time tree
todrop = NULL

# Choose color palette
# Palette gallery: https://pmassicotte.github.io/paletteer_gallery/
color_palette = "pals::kovesi.rainbow_bgyr_35_85_c72"

# Load trees
VR_tree_file = file.path(path_results, paste0(job_name, ".Output.trees"))
vr_trees = read.nexus(VR_tree_file)

# Load functions & map rates on time tree, exports heat mapped tree
source(file.path(path_scripts, "bt_VRPlotting.R"))
vr_mapped = mapped_rates(tree, vr_trees, todrop=todrop)
HeatMap = VR_plot(vr_mapped$mcc_tree, vr_mapped$mapped_data, color_palette)
