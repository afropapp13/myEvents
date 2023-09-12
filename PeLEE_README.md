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

root -b PeLEE_Create1DPlotsTHStack_TopologicalBreakDown.cpp
root -b PeLEE_Create1DPlotsTHStack_InteractionBreakDown.cpp

###########################################################################################################################################################################

