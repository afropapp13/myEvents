. ../myClasses/Constants.sh

declare -a arrRun=("Combined")

declare -a arrCuts=("_NoCuts_PID_NuScore")

# Loop over the run numbers

for RunNumber in "${arrRun[@]}"
do

	AppliedCuts=""

	# Loop over the selection cuts

	for Cut in "${arrCuts[@]}"
	do

		AppliedCuts=${AppliedCuts}${Cut}

		scp $UserID@$UBgpvm:$CCQEPlotPath/$AppliedCuts/TopologicalBreakDown/*THStack_BreakDown_*_${RunNumber}_${UBCode}${AppliedCuts}.pdf \
			./myCCQEPlots/pdf/1D/${UBCode}/${AppliedCuts}/TopologicalBreakDown
			
		scp $UserID@$UBgpvm:$CCQEPlotPath/$AppliedCuts/InteractionBreakDown/*THStack_BreakDown_*_${RunNumber}_${UBCode}${AppliedCuts}.pdf \
			./myCCQEPlots/pdf/1D/${UBCode}/${AppliedCuts}/InteractionBreakDown			

	done

done

# End of the loop over the run numbers
