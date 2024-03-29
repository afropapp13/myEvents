###########################################################################################################################################################################

# switch between series of cuts

root -b PeLEE_script_EventSelection_CV.C 

###########################################################################################################################################################################

# Detector variations
root -b PeLEE_script_EventSelection_Detector_Systematics.C

# GEANT4 variations
root -b PeLEE_script_EventSelection_G4_Systematics.C

# GENIE variations
root -b PeLEE_script_EventSelection_Genie_Systematics.C

# Flux variations
root -b PeLEE_script_EventSelection_Flux_Systematics.C

# MC_Stat variations
root -b PeLEE_script_EventSelection_MC_Stat_Systematics.C

# NuWro Fake Data
root -b PeLEE_script_EventSelection_FakeData.C

###########################################################################################################################################################################

#GENIE v2 Fake Data
root -b PeLEE_script_EventSelection_GENIEv2.C

###########################################################################################################################################################################

# Untuned MC & Twice MEC

root -b PeLEE_script_EventSelection_CV_FakeData.C

root -b PeLEE_script_EventSelection_Detector_Systematics_FakeData.C

root -b PeLEE_script_EventSelection_G4_Systematics_FakeData.C

root -b PeLEE_script_EventSelection_Genie_Systematics_FakeData.C

root -b PeLEE_script_EventSelection_Flux_Systematics_FakeData.C

root -b PeLEE_script_EventSelection_MC_Stat_Systematics_FakeData.C

###########################################################################################################################################################################

# Detailed xsec systematics

#root -b PeLEE_script_EventSelection_DetailedGenie_Systematics.C

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
root -b PRD_InteractionBreakDown.cpp
root -b PRD_SignalBkg_InteractionBreakDown.cpp
#root -b PRL_SuppMat_PID.cpp
# LLR PRD Fig131.pdf
root -b PRD_LLR_ParticleBreakDown.cpp 
# Muon / Proton LLR PRD
root -b PRD_PID_NoCuts.cpp

root -b
.L PeLEE_Create1DPlotsTHStack_TopologicalBreakDown.cpp
PeLEE_Create1DPlotsTHStack_TopologicalBreakDown("Overlay9NuWro")
PeLEE_Create1DPlotsTHStack_TopologicalBreakDown("NoTuneOverlay9")
PeLEE_Create1DPlotsTHStack_TopologicalBreakDown("TwiceMECOverlay9")
PeLEE_Create1DPlotsTHStack_TopologicalBreakDown("GENIEv2Overlay9")

root -b
.L PeLEE_Create1DPlotsTHStack_InteractionBreakDown.cpp
PeLEE_Create1DPlotsTHStack_InteractionBreakDown("Overlay9NuWro")
PeLEE_Create1DPlotsTHStack_InteractionBreakDown("NoTuneOverlay9")
PeLEE_Create1DPlotsTHStack_InteractionBreakDown("TwiceMECOverlay9")
PeLEE_Create1DPlotsTHStack_InteractionBreakDown("GENIEv2Overlay9")

root -b 
.L PeLEE_PrintLatexTables.cpp
PeLEE_PrintLatexTables("",true)
PeLEE_PrintLatexTables("",false,true)
PeLEE_PrintLatexTables("",false,false,true)
PeLEE_PrintLatexTables("",false,false,false,true)
PeLEE_PrintLatexTables("",false,false,false,false,true)
PeLEE_PrintLatexTables("",false,false,false,false,false,true)

###########################################################################################################################################################################

# locally
./PeLEE_DownloadEventRatePlots.sh

###########################################################################################################################################################################

