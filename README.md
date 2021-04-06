# myEvents

# Nominal/CV predictions 
# To be done iteractively on the uboone gpvm's
# The preselected files are located under "/pnfs/uboone/persistent/users/apapadop/mySamples/"+UBCodeVersion+"/PreSelection_"+fWhichSample+"_"+UBCodeVersion+".root";
# The POT files are located under "/pnfs/uboone/persistent/users/apapadop/mySamples/"+UBCodeVersion+"/PreSelection_"+fWhichSample+"_"+UBCodeVersion+"_POT.root"
# The output files will be placed under "/uboone/data/users/apapadop/myEvents/OutputFiles/"+UBCodeVersion+"/"+Cuts+"/STVStudies_"+fWhichSample+Extension+Cuts+".root"

###########################################################################################################################################################################

# switch between series of cuts

root -l script_EventSelection_CV.C 

# switch between PID info / topological breakdown / interaction breakdown

root -l ProduceValidationPlots.cpp

./DownloadEventRatePlots.sh

###########################################################################################################################################################################

# Purity & Efficiency Studies

cd PurityEfficiencyStudies
root -l script_PurityEfficiencyStudies.C
root -l DetermineMaxPurityEfficiency1D.cpp

# Detector variations

root -l script_EventSelection_Detector_Systematics.C

# GEANT4 variations

root -l script_EventSelection_G4_Systematics.C

# GENIE variations

root -l script_EventSelection_Genie_Systematics.C

# Flux variations

root -l script_EventSelection_Flux_Systematics.C

###########################################################################################################################################################################

# Plotting code in ProduceValidationPlots.cpp locally

# 3 plane log likelihood particle breakdown 
gROOT->ProcessLine(".L Chi2PID_BreakDown.cpp++"); gROOT->ProcessLine("Chi2PID_BreakDown()");

# Topological break down
gROOT->ProcessLine(".L Create1DPlotsTHStack_TopologicalBreakDown.cpp++"); gROOT->ProcessLine("Create1DPlotsTHStack_TopologicalBreakDown()");

# Interaction break down
gROOT->ProcessLine(".L Create1DPlotsTHStack_InteractionBreakDown.cpp++"); gROOT->ProcessLine("Create1DPlotsTHStack_InteractionBreakDown()");
