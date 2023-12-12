###########################################################################################################################################################################

# switch between series of cuts

root -b event_selection_cv.cxx 

###########################################################################################################################################################################

# Detector variations
root -b event_selection_det.cxx 

# GEANT4 variations
root -b event_selection_g4.cxx

# GENIE variations
root -b event_selection_xsec.cxx 

# Flux variations
root -b event_selection_flux.cxx 

# MC_Stat variations
root -b event_selection_mc_stat.cxx 

# NuWro Fake Data
root -b event_selection_nuwro_fds.cxx

###########################################################################################################################################################################

# Untuned MC & Twice MEC

root -b event_selection_cv_fds.cxx

root -b event_selection_det_fds.cxx

root -b event_selection_g4_fds.cxx

root -b event_selection_xsec_fds.cxx

root -b event_selection_flux_fds.cxx

root -b event_selection_mc_stat_fds.cxx

###########################################################################################################################################################################

root -b 
.L print_latex_tables.cxx
print_latex_tables("",true)
print_latex_tables("",false,true)
print_latex_tables("",false,false,true)
print_latex_tables("",false,false,false,true)
print_latex_tables("",false,false,false,false,true)

###########################################################################################################################################################################

root -b topological_breakdown.cxx
root -b interaction_breakdown.cxx

###########################################################################################################################################################################

