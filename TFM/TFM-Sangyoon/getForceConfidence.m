function [] = getForceConfidence(pathForMovieData)
% this function calculates force confidence given force parameters and
% output related to movieData
    disp('Loading BEMParams ...')
    tic
    load([pathForMovieData '/TFMPackage/forceField/BEMParams.mat']) % This will load M
    forceNodeMaxima = max(M);
    forceConfidence.pos = forceMesh.p;
    forceConfidence.vec = reshape(forceNodeMaxima,[],2);
    maxCfd = max(forceNodeMaxima);
    forceConfidence.vec = forceConfidence.vec/maxCfd;
    toc

    load([pathForMovieData '/TFMPackage/correctedDisplacementField/displField.mat'])
% figure, imshow(fCfdMapIn{1},[])
% colormap jet
% hold on
% plot(displField(1).pos(:,1),displField(1).pos(:,2),'ko')
    load([pathForMovieData '/movieData.mat']) %this will load MD
    movieData = MD;
    iSDCProc =movieData.getProcessIndex('StageDriftCorrectionProcess',1,1);
    SDCProc=movieData.processes_{iSDCProc};
    pDisp.ChannelIndex = 1;
    refFrame = double(imread(SDCProc.outFilePaths_{2,pDisp.ChannelIndex}));
    firstMask = false(size(refFrame));
    cur_fCfdMap = zeros(size(firstMask));
    ii=1;
    disp('Generating heatmap for force confidence...')
    tic
    [fCfdMapIn, fCmax, fCmin,cropInfo] = generateHeatmapShifted(forceConfidence,displField,0);
    fCfdMapIn{ii} = fCfdMapIn{ii}/max(fCfdMapIn{ii}(:));
    cur_fCfdMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = fCfdMapIn{ii};
    toc
% figure, imshow(cur_fCfdMap,[]), hold on
% plot(displField(1).pos(:,1),displField(1).pos(:,2),'ko')
    disp('Loading existing tractionMaps.mat ...')
    tic
    load([pathForMovieData '/TFMPackage/forceField/tractionMaps.mat'])
    p.OutputDirectory=[pathForMovieData '/TFMPackage/forceField'];
    outputFile{2,1} = [p.OutputDirectory filesep 'tractionMaps.mat'];
    toc
%     fCfdMap = cell(1,1); %force confidence
    fCfdMap = cur_fCfdMap;
    disp('Saving fCfdMaps additionally in tractionMaps.mat ...')
    tic
    save(outputFile{2},'tMap','tMapX','tMapY','fCfdMap'); % need to be updated for faster loading. SH 20141106
    toc
end