. ../myClasses/Constants.sh

declare -a arrRun=("Run1" "Run3")

declare -a arrCuts=("_NoCuts_PID_NuScore")

# Loop over the run numbers

for RunNumber in "${arrRun[@]}"
do

	AppliedCuts=""

	# Loop over the selection cuts

	for Cut in "${arrCuts[@]}"
	do

		AppliedCuts=${AppliedCuts}${Cut}

		scp $UserID@$UBgpvm:$PlotPath/$AppliedCuts/TopologicalBreakDown/THStack_Syst_BreakDown_*_${RunNumber}_${UBCode}${AppliedCuts}.pdf \
			./myPlots/pdf/1D/${UBCode}/${AppliedCuts}/TopologicalBreakDown		

	done

done

# End of the loop over the run numbers
