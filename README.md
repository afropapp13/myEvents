# myEvents

# Nominal predictions

root -l script.C

# Detector variations

root -l script_Detector_Systematics.C



# Plotting code in ProduceValidationPlots.cpp

# 3 plane log likelihood particle breakdown 
gROOT->ProcessLine(".L Chi2PID_BreakDown.cpp++"); gROOT->ProcessLine("Chi2PID_BreakDown()");

# Topological break down
gROOT->ProcessLine(".L Create1DPlotsTHStack_TopologicalBreakDown.cpp++"); gROOT->ProcessLine("Create1DPlotsTHStack_TopologicalBreakDown()");

# Interaction break down
gROOT->ProcessLine(".L Create1DPlotsTHStack_InteractionBreakDown.cpp++"); gROOT->ProcessLine("Create1DPlotsTHStack_InteractionBreakDown()");

# 2D plots # Raw output, not normalized rows to 1 # For that, take a look at mySTVAnalysis/MigrationMatrices.cpp
gROOT->ProcessLine(".L Create2DPlots.cpp++"); gROOT->ProcessLine("Create2DPlots()");
