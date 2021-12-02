# DODCG-sampling-eval
 Scripts for implementing simulation-based evaluation of sampling for Camp Guernsey overlay

# Script file catalog (does not include archived scripts)

00-Sampling_unit_coords.R - Gets spatial coordinates for sampling units to match up with GIS vegetation catagories.
01-Summarize_point_habitats - Tabulates the number of sampling units in each GIS vegetation category.
02-DODCG sample size sim.R - Implements simulations of density estimation scenarios.
02-DODCG trend power sim.R - Implements simulations of trend estimation scenarios.
02-Spp_hab_specializations.R - Categorizes species habitat associations based on relative abundance index.
03-Species_summary.R - Tabulates species summary tables for report.
03-Tabulate spp sigma estimates.R - Tabulates detection distance parameters for classifying species as low vs high detectability.
sample_size_sim_HNdot_identity.jags - JAGS model for density estimation simulations.
trend_sim_HNdot_noRE.jags - JAGS model for trend estimation simulations.