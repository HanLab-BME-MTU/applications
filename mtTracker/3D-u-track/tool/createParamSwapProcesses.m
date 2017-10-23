function processList=createParamSwapProcesses(motherProcess,funParamSwipe)
% Create multiple process copy with different parameter specified by funParamSwipe.
% Each process is equipped with proper outputDirectory

paramToSwipe=fields(funParamSwipe);
processList=[];
for fIdx=1:length(paramToSwipe)
	fieldName=paramToSwipe{fIdx};
	paramValues=getfield(funParamSwipe,paramToSwipe(fIdx));
	if(~isnumeric(paramValues))
		processList=[processList createParamSwapProcesses(motherProcess,paramValues)]
	else
		for pIdx=1:length(paramValues)
			pr=copy(motherProcess);
			funParams=pr.funParams_;
			funParams.setfield(fieldName,paramValues(pIdx));
			funParams.OutputDirectory=[funParams.OutputDirectory '_' fieldName '_' num2str(paramValues)];
			pr.setPara(funParams);
			processList=[processList pr];
			LessParamSwipe=funParamSwipe;
			LessParamSwipe{1:fIdx}=[];
			processList=[processList buildParamSpace(pr,LessParamSwipe)];
		end
	end
end


