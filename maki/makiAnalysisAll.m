analysisStruct = makiSisterConnectionAnalysis('DANUSER',analysisStruct,1,1,1,0);
analysisStruct = makiSisterMotionCoupling('DANUSER',analysisStruct,1,1,1,0);
analysisStruct = makiSepDispSpaceTimeAnalysis('DANUSER',analysisStruct,1,1,1,0);
analysisStruct = makiMetaPlateAnalysis('DANUSER',analysisStruct,1);
% analysisStruct.sepDispSpaceTime = [];
% analysisStruct.metaPlate = [];
