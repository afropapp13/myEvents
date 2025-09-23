###########################################################################################################################################################################

# run just the reco combined files 

root -b mcc9_10_event_selection_cv_reco.cxx 

###########################################################################################################################################################################

# run both reco and truth for all files
./mcc9_10_run_event_selection.sh 

###########################################################################################################################################################################

root -b 
.L mcc9_10_print_latex_tables.cxx
# print stats
mcc9_10_print_latex_tables("",true)
# print purity / efficiency
mcc9_10_print_latex_tables("",false,true)
# cosmic/dirt contamination
mcc9_10_print_latex_tables("",false,false,true)
# interaction breakdown
mcc9_10_print_latex_tables("",false,false,false,true)

###########################################################################################################################################################################
