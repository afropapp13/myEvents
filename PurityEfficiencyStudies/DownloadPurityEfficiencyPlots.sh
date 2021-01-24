. ../../myClasses/Constants.sh

export OutPutDir=/uboone/data/users/$UserID/mySTVAnalysis/myPlots/$UBCode/

scp $UserID@$UBgpvm:$OutPutDir/*Purity*.pdf ../../myPreSelection/myPlots/$UBCode/
scp $UserID@$UBgpvm:$OutPutDir/*Efficiency*.pdf ../../myPreSelection/myPlots/$UBCode/
scp $UserID@$UBgpvm:$OutPutDir/*Product*.pdf ../../myPreSelection/myPlots/$UBCode/
