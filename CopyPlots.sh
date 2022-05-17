. ../myClasses/Constants.sh

export InPutDirNoCuts=myPlots/pdf/1D/v08_00_00_52/_NoCuts/TopologicalBreakDown/
export InPutDir=myPlots/pdf/1D/v08_00_00_52/_NoCuts_PID_NuScore/InteractionBreakDown/
export OutPutDir=/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_KinematicImbalance/Figures
export OutPutDirPublicNote=/home/afroditi/Dropbox/Apps/Overleaf/MicroBooNE_Neutrino2022_PublicNote/Figures

declare -a arrPlots=(
"THStack_BreakDown_RecoDeltaPTPlot_Combined_"${UBCode}"_NoCuts_PID_NuScore.pdf"
"THStack_BreakDown_RecoDeltaAlphaTPlot_Combined_"${UBCode}"_NoCuts_PID_NuScore.pdf"
"THStack_BreakDown_RecoDeltaPtxPlot_Combined_"${UBCode}"_NoCuts_PID_NuScore.pdf"
"THStack_BreakDown_RecoDeltaPtyPlot_Combined_"${UBCode}"_NoCuts_PID_NuScore.pdf"
)

for plot in "${arrPlots[@]}"
do

	##############################################################################

	cp ${InPutDir}/${plot}	${OutPutDir}
	cp ${InPutDir}/${plot}	${OutPutDirPublicNote}	
	
	##############################################################################

done

cp ${InPutDirNoCuts}/../PID_BreakDown_RecoLLRPIDPlot_NoCuts_Combined_v08_00_00_52.pdf		${OutPutDir}
cp ${InPutDirNoCuts}/THStack_BreakDown_RecoProtonLLRPIDPlot_Combined_v08_00_00_52_NoCuts.pdf	${OutPutDir}
