# Public simulation design inputs for the SMART_HM simulations.
#
# These values define the data-generating model used by the public simulation
# scripts. They are intentionally stored as code so the repository does not
# depend on private application data or generated calibration result folders.

simulation_design_params <- conditional_design_defaults()

# Reference baseline outcome value used when translating interpretable design
# inputs into conditional hurdle GLM coefficients.
simulation_Y0_ref <- 12
