function nParam = nParamModel(modelType,fitMethod)

maxType = max(modelType);

switch(fitMethod)
    
    case 1 % std. 3D
        nParamBez = [3;((1:maxType)'+1)*2+2];
        
    case 2 % std. 2D
        nParamBez = [2;((1:maxType)'+1)+2];
        
    case 3 % snakes 3D
        nParamBez = ((0:maxType)'+1)*3;
        
    case 4 % snakes 2D
        nParamBez = ((0:maxType)'+1)*2;

end

nParamVar = 1;
nParam = nParamBez+nParamVar;
nParam = nParam(modelType+1);

end

