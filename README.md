###########################################################################################################################################################################

# run just the reco combined files 

root -b event_selection_cv_reco.cxx 

###########################################################################################################################################################################

# run both reco and truth for all files
./run_event_selection.sh 

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

# NuWro
root -b
.L topological_breakdown.cxx 
topological_breakdown("Overlay9NuWro")

root -b
.L interaction_breakdown.cxx 
interaction_breakdown("Overlay9NuWro")


###########################################################################################################################################################################

