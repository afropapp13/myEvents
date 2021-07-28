. ../../myClasses/Constants.sh

export OutPutDir=/uboone/data/users/$UserID/mySTVAnalysis/myCCQEPlots/$UBCode/

scp $UserID@$UBgpvm:$OutPutDir/*Purity*.pdf ../../myCCQEPreSelection/myCCQEPlots/$UBCode/
scp $UserID@$UBgpvm:$OutPutDir/*Efficiency*.pdf ../../myCCQEPreSelection/myCCQEPlots/$UBCode/
scp $UserID@$UBgpvm:$OutPutDir/*Product*.pdf ../../myCCQEPreSelection/myCCQEPlots/$UBCode/
