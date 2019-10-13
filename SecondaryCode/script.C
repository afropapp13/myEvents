{
gROOT->ProcessLine(".L ../SecondaryCode/TruncMean.cxx");
gROOT->ProcessLine(".L anatree.C");
gROOT->ProcessLine("anatree().Loop()");
//gROOT->ProcessLine(".q");
};
