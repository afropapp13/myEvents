. ../myClasses/Constants.sh

export NoCuts=_NoCuts
#export Cuts=_NoCuts_NuScore_ThreePlaneLogChi2_Collinearity
export Cuts=_NoCuts_PID_NuScore

export OutPutDir=/uboone/data/users/$UserID/myEvents/OutputFiles/$UBCode

declare -a arrRun=("Run1" "Run3")
#declare -a arrRun=("Run3")
declare -a arrCuts=("_NoCuts" "_PID" "_NuScore")

# Loop over the run numbers

for RunNumber in "${arrRun[@]}"
do

	##############################################################################

	# Truth level info

	scp $UserID@$UBgpvm:$OutPutDir/TruthSTVAnalysis_Overlay9_${RunNumber}_${UBCode}.root ./OutputFiles/$UBCode/TruthSTVAnalysis_Overlay9_${RunNumber}_${UBCode}.root

	##############################################################################

	AppliedCuts=""

	# Loop over the selection cuts

	for Cut in "${arrCuts[@]}"
	do

		AppliedCuts=${AppliedCuts}${Cut}

		# Reco level info # All cuts applied

		scp $UserID@$UBgpvm:$OutPutDir/$AppliedCuts/STVStudies_Overlay9_$RunNumber$AppliedCuts.root \
			./OutputFiles/$UBCode/$AppliedCuts/STVStudies_Overlay9_$RunNumber$AppliedCuts.root

		scp $UserID@$UBgpvm:$OutPutDir/$AppliedCuts/STVStudies_OverlayDirt9_$RunNumber$AppliedCuts.root \
			./OutputFiles/$UBCode/$AppliedCuts/STVStudies_OverlayDirt9_$RunNumber$AppliedCuts.root

		scp $UserID@$UBgpvm:$OutPutDir/$AppliedCuts/STVStudies_BeamOn9_$RunNumber$AppliedCuts.root \
			./OutputFiles/$UBCode/$AppliedCuts/STVStudies_BeamOn9_$RunNumber$AppliedCuts.root

		scp $UserID@$UBgpvm:$OutPutDir/$AppliedCuts/STVStudies_ExtBNB9_$RunNumber$AppliedCuts.root \
			./OutputFiles/$UBCode/$AppliedCuts/STVStudies_ExtBNB9_$RunNumber$AppliedCuts.root

	done

done

# End of the loop over the run numbers
