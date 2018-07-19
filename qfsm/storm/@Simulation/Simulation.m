classdef Simulation < handle
    
    % ---------------------
    % Class providing various methods to generate synthetic data.
    % Pascal Bérard, February 2012
    % ---------------------
    
    properties (GetAccess = 'private',SetAccess = 'private')
        samplingType = 'regular'; % 'random';
        data;
    end
    
    methods
        % Constructor
        function obj = Simulation(data)
            if nargin>0
                obj.data = data;
            end
        end
        
        addGaussianNoise(obj,sigmaNoise);
        addGaussianNoiseGradientX(obj,sigmaNoiseStart,sigmaNoiseEnd,xStart,xEnd);
        addRandomNoise(obj,nNoisePoints);
        t = bezier(obj,controlPoints,nSamples);
        t = bezierWithDensity(obj,controlPoints,density);
        t = bezierWithOrthogonalGaussianNoise2D(obj,controlPoints,nSamples,sigmaNoise);
        filament(obj,controlPoints,nSamples,sigmaNoise,sourceSpacing,fractionOfDamagedSources);
        t = filamentFromPDF(obj,controlPoints,sigmaNoise,pdf);
        setDomain(obj,roiPosition,roiSize);
        line(obj,startPoint,endPoint,nSamples);
        lineWithDensity(obj,startPoint,endPoint,density);
        lineWithDensityGradient(obj,startPoint,endPoint,densityStart,densityEnd);
        meshDensity(obj,height,length,nHorzLines,squeezeFactor,density,sigmaNoise);
        meshNoise(obj,height,length,nHorzLines,density,sigmaNoiseStart,sigmaNoiseEnd);
        meshSampling(obj,height,length,nHorzLines,densityStart,densityEnd,sigmaNoise);
        modelToPoints(obj,model,sigmaNoise,density);
        modelToPointsFromPDF(obj,model,sigmaNoise,pdf);
        setSamplingToRandom(obj);
        setSamplingToRegular(obj);
        
    end
    
end

