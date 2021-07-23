# myEvents

# Nominal/CV predictions 
# To be done iteractively on the uboone gpvm's
# The preselected files are located under "/pnfs/uboone/persistent/users/apapadop/mySamples/"+UBCodeVersion+"/PreSelection_"+fWhichSample+"_"+UBCodeVersion+".root";
# The POT files are located under "/pnfs/uboone/persistent/users/apapadop/mySamples/"+UBCodeVersion+"/PreSelection_"+fWhichSample+"_"+UBCodeVersion+"_POT.root"
# The output files will be placed under "/uboone/data/users/apapadop/myEvents/OutputFiles/"+UBCodeVersion+"/"+Cuts+"/STVStudies_"+fWhichSample+Extension+Cuts+".root"

###########################################################################################################################################################################

# switch between series of cuts

root -b PeLEE_script_EventSelection_CV.C 

###########################################################################################################################################################################

# Purity & Efficiency Studies

# on gpvm's
cd PurityEfficiencyStudies
root -b DetermineMaxPurityEfficiency1D.cpp

# (locally)
cd PurityEfficiencyStudies
./DownloadPurityEfficiencyPlots.sh

###########################################################################################################################################################################

root -b PeLEE_Chi2PID_BreakDown.cpp

root -b PeLEE_Create1DPlotsTHStack_TopologicalBreakDown.cpp

root -b PeLEE_Create1DPlotsTHStack_InteractionBreakDown.cpp

root -b PeLEE_PrintLatexTables.cpp

# locally
./PeLEE_DownloadEventRatePlots.sh

###########################################################################################################################################################################

# Detector variations

root -b PeLEE_script_EventSelection_Detector_Systematics.C

# GEANT4 variations

root -b PeLEE_script_EventSelection_G4_Systematics.C

# GENIE variations

root -b PeLEE_script_EventSelection_Genie_Systematics.C

# Flux variations

root -b PeLEE_script_EventSelection_Flux_Systematics.C

###########################################################################################################################################################################

# NuWro fake data (only Run 1 for now)

root -b
.L PeLEE_Create1DPlotsTHStack_TopologicalBreakDown.cpp
PeLEE_Create1DPlotsTHStack_TopologicalBreakDown("Overlay9NuWro")

root -b
.L PeLEE_Create1DPlotsTHStack_InteractionBreakDown.cpp
PeLEE_Create1DPlotsTHStack_InteractionBreakDown("Overlay9NuWro")

# NoTune GENIE Fake Data

root -b PeLEE_script_EventSelection_CV_FakeData.C 

root -b
.L PeLEE_Create1DPlotsTHStack_TopologicalBreakDown.cpp
PeLEE_Create1DPlotsTHStack_TopologicalBreakDown("NoTuneOverlay9")

root -b
.L Create1DPlotsTHStack_InteractionBreakDown.cpp
Create1DPlotsTHStack_InteractionBreakDown("NoTuneOverlay9")

