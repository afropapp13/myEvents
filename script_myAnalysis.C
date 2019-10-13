{
gROOT->ProcessLine(".L ../../MyClasses/Tools.cxx+");
gROOT->ProcessLine(".L myAnalysis.C+");
gROOT->ProcessLine("myAnalysis().Loop()");
gROOT->ProcessLine(".q");
};
