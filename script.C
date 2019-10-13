{
gROOT->ProcessLine(".L t.C+");
gROOT->ProcessLine("t().Loop()");
gROOT->ProcessLine(".q");
};
