{
gROOT->ProcessLine(".L ../../MyClasses/Tools.cxx+");
gROOT->ProcessLine(".L myTrueAnalysis.C+");
gROOT->ProcessLine("myTrueAnalysis().Loop()");
gROOT->ProcessLine(".q");
};
