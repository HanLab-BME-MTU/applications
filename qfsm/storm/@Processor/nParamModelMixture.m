function nParam = nParamModelMixture(modelType1,modelType2,fitMethod)

nParamModels = Processor.nParamModel(modelType1,fitMethod)+Processor.nParamModel(modelType2,fitMethod);
nParamMix = 2; % The mixture parameters
nParam = nParamModels+nParamMix;

end
