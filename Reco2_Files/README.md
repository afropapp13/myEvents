# myEvents

# Nominal/CV predictions 
# To be done iteractively on the uboone gpvm's
# The preselected files are located under "/pnfs/uboone/persistent/users/apapadop/mySamples/"+UBCodeVersion+"/PreSelection_"+fWhichSample+"_"+UBCodeVersion+".root";
# The POT files are located under "/pnfs/uboone/persistent/users/apapadop/mySamples/"+UBCodeVersion+"/PreSelection_"+fWhichSample+"_"+UBCodeVersion+"_POT.root"
# The output files will be placed under "/uboone/data/users/apapadop/myEvents/OutputFiles/"+UBCodeVersion+"/"+Cuts+"/STVStudies_"+fWhichSample+Extension+Cuts+".root"

###########################################################################################################################################################################

# switch between series of cuts

root -l script_EventSelection_CV.C 

###########################################################################################################################################################################

# Purity & Efficiency Studies

# on gpvm's
cd PurityEfficiencyStudies
root -l DetermineMaxPurityEfficiency1D.cpp

# (locally)
cd PurityEfficiencyStudies
./DownloadPurityEfficiencyPlots.sh

###########################################################################################################################################################################

root -l Chi2PID_BreakDown.cpp

root -l Create1DPlotsTHStack_TopologicalBreakDown.cpp

root -l Create1DPlotsTHStack_InteractionBreakDown.cpp

root -l PrintLatexTables.cpp

./DownloadEventRatePlots.sh

###########################################################################################################################################################################

# Detector variations

root -l script_EventSelection_Detector_Systematics.C

# GEANT4 variations

root -l script_EventSelection_G4_Systematics.C

# GENIE variations

root -l script_EventSelection_Genie_Systematics.C

# Flux variations

root -l script_EventSelection_Flux_Systematics.C

###########################################################################################################################################################################

# NuWro fake data (only Run 1 for now)

root -l
.L Create1DPlotsTHStack_TopologicalBreakDown.cpp
Create1DPlotsTHStack_TopologicalBreakDown("Overlay9NuWro")

root -l
.L Create1DPlotsTHStack_InteractionBreakDown.cpp
Create1DPlotsTHStack_InteractionBreakDown("Overlay9NuWro")

# NoTune GENIE Fake Data

root -l script_EventSelection_CV_FakeData.C 

root -l
.L Create1DPlotsTHStack_TopologicalBreakDown.cpp
Create1DPlotsTHStack_TopologicalBreakDown("NoTuneOverlay9")

root -l
.L Create1DPlotsTHStack_InteractionBreakDown.cpp
Create1DPlotsTHStack_InteractionBreakDown("NoTuneOverlay9")

