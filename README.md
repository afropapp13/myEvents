###########################################################################################################################################################################

# run just the reco combined files 

root -b event_selection_cv_reco.cxx 

###########################################################################################################################################################################

# run both reco and truth for all files
./run_event_selection.sh 

###########################################################################################################################################################################

root -b 
.L print_latex_tables.cxx
# print stats
print_latex_tables("",true)
# print purity / efficiency
print_latex_tables("",false,true)
# cosmic/dirt contamination
print_latex_tables("",false,false,true)
# interaction breakdown
print_latex_tables("",false,false,false,true)

###########################################################################################################################################################################

root -b print_1d_slices.cxx
root -b print_2d_slices.cxx
root -b fluxes.cxx
root -b print_theta_vis_mean_std.cxx
root -b dune_projection.cxx
root -b data_dune_projection.cxx

###########################################################################################################################################################################

