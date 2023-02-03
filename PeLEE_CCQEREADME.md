# myEvents

# Nominal/CV predictions 
# To be done iteractively on the uboone gpvm's
# The preselected files are located under "/pnfs/uboone/persistent/users/apapadop/mySamples/"+UBCodeVersion+"/PreSelection_"+fWhichSample+"_"+UBCodeVersion+".root";
# The POT files are located under "/pnfs/uboone/persistent/users/apapadop/mySamples/"+UBCodeVersion+"/PreSelection_"+fWhichSample+"_"+UBCodeVersion+"_POT.root"
# The output files will be placed under "/uboone/data/users/apapadop/myEvents/OutputFiles/"+UBCodeVersion+"/"+Cuts+"/STVStudies_"+fWhichSample+Extension+Cuts+".root"

###########################################################################################################################################################################

root -b PeLEE_script_CCQEEventSelection_CV.C 

root -b PeLEE_CCQECreate1DPlotsTHStack_InteractionBreakDown.cpp

root -b PeLEE_CCQECreate1DPlotsTHStack_TopologicalBreakDown.cpp

# locally
./PeLEE_CCQEDownloadEventRatePlots.sh

###########################################################################################################################################################################

# Detector variations

root -b PeLEE_script_CCQEEventSelection_Detector_Systematics.C

# GEANT4 variations

root -b PeLEE_script_CCQEEventSelection_G4_Systematics.C

# GENIE variations

root -b PeLEE_script_CCQEEventSelection_Genie_Systematics.C

# Flux variations

root -b PeLEE_script_CCQEEventSelection_Flux_Systematics.C

# MC_Stat variations
root -b PeLEE_script_CCQEEventSelection_MC_Stat_Systematics.C

#NuWro
root -b PeLEE_script_CCQEEventSelection_FakeData.C
