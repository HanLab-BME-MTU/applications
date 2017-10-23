function [maxSpeed,processOverlay]=trackabilityAnalysis(MD,processDetection,processProj,processBuildRef)  

% TODO: Generic Package
    MD.getPackage(160).displayProcessInfo
    processProj=MD.getPackage(160).getProcess(8);
    processDetection=MD.getPackage(160).getProcess(2);
    detections=load(processDetection.outFilePaths_{1});
    detections=Detections(detections.movieInfo);
    [maxSpeed,~,densities]=estimateTrackability(detections,1/MD.timeInterval_,'debugMode','AmbiguityHeuristic');
    colorIndex=cellfun(@(m) ceil(254*mat2gray(m,[0.5,0.8]))+1,maxSpeed,'unif',0);
    %colorIndex=cellfun(@(m) ceil(254*mat2gray(m,[0,10]))+1,densities,'unif',0);
    processBuildRef=MD.getPackage(160).getProcess(4);
    refs=load(processBuildRef.outFilePaths_{1}); refs=refs.refs;
    ROIs=load(processBuildRef.outFilePaths_{2}); ROIs=ROIs.ROI;
    myColormap=255*flipud(jet(256));
    processOverlay=ExternalProcess(MD,'Overlay');
    overlayProjDetectionMovie(processProj,'detections',  refs(1,2).applyBase(detections,''), ...
        'colorIndx',colorIndex, ...
        'process',processOverlay, ...
        'colormap',myColormap,'name','trackability','processFrames',[1:5]);
    MD.addPackage(GenericPackage({processProj,processOverlay},[],'name_','trackabilityAnalysis'));