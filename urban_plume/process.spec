
process env                     # environment state
suffix env                      # filename suffix for output

process gas                     # gas state
suffix gas                      # filename suffix for output

process aero                    # binned aerosol state
suffix aero                     # filename suffix for output

process kappa                   # kappa-derived critical supersaturation
suffix kappa_crit_ss            # filename suffix for output
n_step 100                      # number of steps for histogram
min 1e-5                        # min supersaturation (1)
max 5e-2                        # max supersaturation (1)

process comp                    # scalar composition
suffix comp_bc                  # filename suffix for output
n_step 100                      # number of steps for histogram
min 0                           # minimum ratio (1)
max 1                           # maximum ratio (1)
a_species SO4_a NO3_a Cl_a NH4_a MSA_a ARO1_a ARO2_a ALK1_a OLE1_a API1_a API2_a LIM1_a LIM2_a CO3_a Na_a Ca_a OIN_a OC_a # first species list
b_species BC_a                  # second species list

#process comp                    # scalar composition
#suffix comp_nh4                 # filename suffix for output
#n_step 100                      # number of steps for histogram
#min 0                           # minimum ratio (1)
#max 1                           # maximum ratio (1)
#a_species SO4_a NO3_a BC_a      # first species list
#b_species NH4_a                 # second species list

#process n_orig_part             # number of original particles
#suffix n_orig_part              # filename suffix for output
#min 0                           # min original particles
#max 100                         # max original particles

process optic_absorb            # optical absorption
suffix optic_absorb             # filename suffix for output
n_step 100                      # number of steps for histogram
min 0                           # min absorption cross-section (m^2)
max 1e-5                        # max absorption cross-section (m^2)

process optic_scatter           # optical scattering
suffix scatter                  # filename suffix for output
n_step 100                      # number of steps for histogram
min 0                           # min scattering cross-section (m^2)
max 1e-5                        # max scattering cross-section (m^2)

process optic_extinct           # optical extinction
suffix extinct                  # filename suffix for output
n_step 100                      # number of steps for histogram
min 0                           # min extinction cross-section (m^2)
max 1e-5                        # max extinction cross-section (m^2)
