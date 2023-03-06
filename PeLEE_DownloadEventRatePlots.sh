. ../myClasses/Constants.sh

declare -a arrRun=("Run1" "Run2" "Run3" "Run4a" "Run4b" "Run4aRutgers" "Combined")

declare -a arrCuts=("_NoCuts" "_PID_NuScore")

#Run4a validation plots
scp $UserID@$UBgpvm:${PlotPath}/*.pdf ./myPlots/pdf/1D/${UBCode}/

#PID plots
scp $UserID@$UBgpvm:${PlotPath}/_NoCuts/*.pdf ./myPlots/pdf/1D/${UBCode}/_NoCuts/
scp apapadop@$UBgpvm.fnal.gov:/uboone/data/users/apapadop/mySTVAnalysis/myPlots/v08_00_00_52/PRL_SuppMat_RecoProtonLLRPIDPlot_Combined_v08_00_00_52_NoCuts.pdf ./myPlots/pdf/1D/v08_00_00_52/

# Loop over the run numbers

for RunNumber in "${arrRun[@]}"
do

	AppliedCuts=""

	# Loop over the selection cuts

	for Cut in "${arrCuts[@]}"
	do

		AppliedCuts=${AppliedCuts}${Cut}

		scp $UserID@$UBgpvm:$PlotPath/$AppliedCuts/TopologicalBreakDown/*THStack_BreakDown_*_${RunNumber}_${UBCode}${AppliedCuts}.pdf \
			./myPlots/pdf/1D/${UBCode}/${AppliedCuts}/TopologicalBreakDown
			
		scp $UserID@$UBgpvm:$PlotPath/$AppliedCuts/InteractionBreakDown/*THStack_BreakDown_*_${RunNumber}_${UBCode}${AppliedCuts}.pdf \
			./myPlots/pdf/1D/${UBCode}/${AppliedCuts}/InteractionBreakDown
			
		scp $UserID@$UBgpvm:$PlotPath/$AppliedCuts/InteractionBreakDown/PRD_Reco*_${RunNumber}_${UBCode}${AppliedCuts}.pdf \
			./myPlots/pdf/1D/${UBCode}/${AppliedCuts}/InteractionBreakDown						

	done

done

# End of the loop over the run numbers
