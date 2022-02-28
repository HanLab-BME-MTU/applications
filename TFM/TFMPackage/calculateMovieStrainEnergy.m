function []=calculateMovieStrainEnergy(movieData, varargin)
%function [SE_FOV,SE_FAs,SE_CellSeg]=quantifyMovieStrainEnergy(curMD)
%quantifies from traction map the strain energy for entire field of view
%(FOV), segmented FAs, and cell segmentation when the cell segmentation
%information is there.
% Unit is in femto Joule (1e-15 J)
% Sangyoon Han March, 2016
% calculateMovieStrainEnergy calculate the strain energy and total force
% out of the force field and the displacement field
%
% calculateMovieStrainEnergy 
%
% SYNOPSIS calculateMovieStrainEnergy(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%

% Sangyoon J. Han, May 2017

disp(['Working on ' movieData.getFullPath '...'])
%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous stage drift processes                                                                     
iProc = movieData.getProcessIndex('StrainEnergyCalculationProcess',1,0);
%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(StrainEnergyCalculationProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
strainEnergyCalcProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(strainEnergyCalcProc,paramsIn);

%% Output setup
%% Backup the original vectors to backup folder
if exist(p.OutputDirectory,'dir')
    contents=dir(p.OutputDirectory);
    if any(arrayfun(@(x) ~x.isdir,contents)) % There should be something inside. Otherwise, backing up
        disp('Backing up the original data')
        ii = 1;
        backupFolder = [p.OutputDirectory ' Backup ' num2str(ii)];
        while exist(backupFolder,'dir')
            backupFolder = [p.OutputDirectory ' Backup ' num2str(ii)];
            ii=ii+1;
        end
        mkdir(backupFolder);
        copyfile(p.OutputDirectory, backupFolder,'f')
    end
end
mkClrDir(p.OutputDirectory);
%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',strainEnergyCalcProc.getName());
end
nFrames = movieData.nFrames_;

%% Load the forcefield
TFMPackage = movieData.getPackage(movieData.getPackageIndex('TFMPackage'));
iForceFieldProc = 4;
forceFieldProc=TFMPackage.processes_{iForceFieldProc};

% tractionImgFolder=p.OutputDirectory;
% imgPath = [tractionImgFolder filesep 'imgs'];
% dataPath = [tractionImgFolder filesep 'data'];
% tifForcePath = [imgPath filesep 'tifForce'];
% tifBSensorPath = [imgPath filesep 'tifBSensor'];
% epsPath = [imgPath filesep 'eps'];
% figPath = [imgPath filesep 'figs'];
% if ~exist(tifForcePath,'dir')
%     system(['mkdir -p ' imgPath]);
%     system(['mkdir -p ' dataPath]);
%     system(['mkdir -p ' epsPath]);
%     system(['mkdir -p ' figPath]);
%     system(['mkdir -p ' tifForcePath]);
%     system(['mkdir -p ' tifBSensorPath]);
%     mkdir(imgPath);
%     mkdir(dataPath);
%     mkdir(epsPath);
%     mkdir(figPath);
%     mkdir(tifForcePath);
%     mkdir(tifBSensorPath);
% end
% Set up the output directories
outputFile{1,1} = [p.OutputDirectory filesep 'strainEnergyInFOV.mat'];
outputFile{1,2} = [p.OutputDirectory filesep 'strainEnergyInCell.mat'];
outputFile{1,3} = [p.OutputDirectory filesep 'forceBlobs.mat'];
if p.exportCSV
    outputFile{1,4} = [p.OutputDirectory filesep 'strainEnergyInFOV.csv'];
    outputFile{1,5} = [p.OutputDirectory filesep 'totalAvgForceInFOV.csv'];
    outputFile{1,6} = [p.OutputDirectory filesep 'strainEnergyInCell.csv'];
    outputFile{1,7} = [p.OutputDirectory filesep 'totalForceInCell.csv'];
    outputFile{1,8} = [p.OutputDirectory filesep 'strainEnergyInForceBlobs.csv'];
    outputFile{1,9} = [p.OutputDirectory filesep 'indivForceInForceBlobs.csv'];
end    
outputFile{1,10} = [p.OutputDirectory filesep 'MovieStrainEnergyData_all.mat'];

maskFolder = [p.OutputDirectory filesep 'BandMasks'];
mkdir(maskFolder)
fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
numStr = @(frame) num2str(frame,fString);
maskIndPath = @(frame) [maskFolder filesep 'mask' numStr(frame) '.mat'];

strainEnergyCalcProc.setOutFilePaths(outputFile);

logMsg='Loading traction map...';
if feature('ShowFigureWindows'), waitbar(0,wtBar,sprintf(logMsg)); end

% tMap=forceFieldProc.loadChannelOutput('output','tMap');
% tMap=forceFieldProc.loadChannelOutput('output','tMapUnshifted','iFrame',ii);
% try
%     tMapObj = tractionMaps.tMap; % this is currently in Pa per pixel (1pix x 1pix)
%     fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
%     numStr = @(frame) num2str(frame,fString);
%     outputDir = fullfile(forceFieldProc.funParams_.OutputDirectory,'tractionMaps');
%     mkdir(outputDir);
%     outFileTMap = @(frame) [outputDir filesep 'tractionMap' numStr(frame) '.mat'];
%     for ii=nFrames:-1:1
%         cur_tMapObj = load(outFileTMap(ii),tMapObj.eachTMapName);
%         tMap{ii} = cur_tMapObj.cur_tMap;
%         progressText((nFrames-ii)/nFrames,'Traction map loading') % Update text
%     end
% catch
%     % Set up the output directories
%     outputDir = fullfile(forceFieldProc.funParams_.OutputDirectory,'tractionMaps');
%     fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
%     numStr = @(frame) num2str(frame,fString);
%     outFileTMap=@(frame) [outputDir filesep 'tractionMap' numStr(frame) '.mat'];
%     tMap=cell(1,nFrames);
% end
yModulus = forceFieldProc.funParams_.YoungModulus;
%% Load the displfield
iCorrectedDisplFieldProc = 3;
CorrectedDisplFieldProc=TFMPackage.processes_{iCorrectedDisplFieldProc};
% logMsg='Loading displacement map...';
% if feature('ShowFigureWindows'), waitbar(0,wtBar,sprintf(logMsg)); end
if ~isempty(CorrectedDisplFieldProc)
    dispProc = CorrectedDisplFieldProc;
%     dMap=CorrectedDisplFieldProc.loadChannelOutput('output','dMapUnshifted');
%     try
%         displMaps=load(CorrectedDisplFieldProc.outFilePaths_{2});
%         dMapObj=displMaps.dMap; % this is currently in pix
%         fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
%         numStr = @(frame) num2str(frame,fString);
%         outputDir = fullfile(CorrectedDisplFieldProc.funParams_.OutputDirectory,'displMaps');
%         outFileDMap = @(frame) [outputDir filesep 'displMap' numStr(frame) '.mat'];
%         for ii=nFrames:-1:1
%             cur_dMapObj = load(outFileDMap(ii),dMapObj.eachDMapName);
%             dMap{ii} = cur_dMapObj.cur_dMap;
%             progressText((nFrames-ii)/nFrames,'One-time displacement map loading') % Update text
%         end
%     catch
%         disp('Assuming displacement proportional to tMap because CorrectedDisplFieldProc was empty')
%         dMap=cellfun(@(x) x/(yModulus),tMap,'UniformOutput',false);
%     end
else
    disp('No iCorrectedDisplFieldProc found. Trying to use displacementField from step 2...')
    iCalculatedDisplFieldProc = 2;
    CalculatedDisplFieldProc=TFMPackage.processes_{iCalculatedDisplFieldProc};
    if ~isempty(CalculatedDisplFieldProc)
        dMap=CalculatedDisplFieldProc.loadChannelOutput('output','dMap');
    else
        disp('Assuming displacement proportional to tMap because there were no dmaps found')
        dMap=cellfun(@(x) x/(yModulus),tMap,'UniformOutput',false);
    end
end
%% Calculate strain energy for FOV
% gridSpacing=1;
pixSize_mu=movieData.pixelSize_*1e-3; % in um/pixel
% factorConvert=gridSpacing^2*pixSize_mu^3/10^6;
% factorConvert=gridSpacing^2*(pixSize_mu*1e-6)^3;
areaConvert=pixSize_mu^2; % in um2/pixel
SE_FOV=struct('SE',zeros(nFrames,1),'area',zeros(nFrames,1),'SEDensity',zeros(nFrames,1));
SE_FOV_SE = zeros(nFrames,1);
SE_FOV_area = zeros(nFrames,1);
SE_FOV_SEDensity = zeros(nFrames,1);
totalForceFOV = zeros(nFrames,1);
avgTractionFOV = zeros(nFrames,1);

SE_Cell=struct('SE',zeros(nFrames,1),'area',zeros(nFrames,1),'SEDensity',zeros(nFrames,1),...
    'SE_peri',zeros(nFrames,1),'SE_inside',zeros(nFrames,1),'SEDensityPeri',zeros(nFrames,1),'SEDensityInside',zeros(nFrames,1));
SE_Cell_SE = zeros(nFrames,1); % in femto-Joule=1e15*(N*m)
SE_Cell_area = zeros(nFrames,1); % this is in um2
SE_Cell_SEDensity = zeros(nFrames,1); % J/m2
SE_Cell_SE_peri = zeros(nFrames,1); % in femto-Joule=1e15*(N*m)
SE_Cell_SE_inside = zeros(nFrames,1); % in femto-Joule=1e15*(N*m)
SE_Cell_SEDensityPeri = zeros(nFrames,1); % J/m2
SE_Cell_SEDensityInside = zeros(nFrames,1); % J/m2

totalForceCell = zeros(nFrames,1);
avgTractionCell= zeros(nFrames,1);
totalForceCellPeri= zeros(nFrames,1);
totalForceCellInside= zeros(nFrames,1);
avgTractionCellPeri= zeros(nFrames,1);
avgTractionCellInside= zeros(nFrames,1);

totalDispCell = zeros(nFrames,1);
totalDispCellPeri = zeros(nFrames,1);
totalDispCellInside = zeros(nFrames,1);
avgDispCell = zeros(nFrames,1);
avgDispCellPeri = zeros(nFrames,1);
avgDispCellInside = zeros(nFrames,1);

SE_Blobs=struct('SE',zeros(nFrames,1),'nFA',zeros(nFrames,1),'areaFA',zeros(nFrames,1),...
    'SEDensity',zeros(nFrames,1),'avgFAarea',zeros(nFrames,1),'avgSEperFA',zeros(nFrames,1));
totalForceBlobs = struct('force',zeros(nFrames,1),'avgTraction',[],...
    'maxTraction',[],'forceBlobPixelIdxList',[]);

SE_Blobs_SE = zeros(nFrames,1); % this is in femto-Joule
SE_Blobs_nFA = zeros(nFrames,1);
SE_Blobs_areaFA = zeros(nFrames,1); % in um2
SE_Blobs_avgFAarea = zeros(nFrames,1);
SE_Blobs_avgSEperFA = zeros(nFrames,1); % still in femto-J
SE_Blobs_SEDensity = zeros(nFrames,1); % J/m2

totalForceBlobs_force= zeros(nFrames,1); % in nN %*(pixSize_mu*1e-6)^2*1e9; % in nN
totalForceBlobs_avgTraction= cell(nFrames,1); % in Pa
totalForceBlobs_maxTraction= cell(nFrames,1); % in Pa
totalForceBlobs_forceBlobPixelIdxList= cell(nFrames,1);
totalForceBlobs_avgTractionCell= cell(nFrames,1); % in Pa

%% Get SDC
% Cell Boundary Mask 
iTFMPack = movieData.getPackageIndex('TFMPackage');
TFMPack=movieData.packages_{iTFMPack}; iSDCProc=1;
SDCProc=TFMPack.processes_{iSDCProc};
if ~isempty(SDCProc)
    try
        iBeadChan=SDCProc.funParams_.iBeadChannel;
    catch
        iBeadChan=1;
    end
else
    iBeadChan=1;
end
% iSDCProc =movieData.getProcessIndex('StageDriftCorrectionProcess',1,1); 
existSDC=false;
if ~isempty(SDCProc)
%     SDCProc=movieData.processes_{iSDCProc};
    s = load(SDCProc.outFilePaths_{3,iBeadChan},'T');    
    T = s.T;
    existSDC=true;
end
%% Get the cell boundary from the best mask
iMaskProcess=movieData.getProcessIndex('MaskIntersectionProcess');
existMask=false;
if ~isempty(iMaskProcess)
    maskProc = movieData.getProcess(iMaskProcess);
    existMask = true;
else
    iMaskProcess=movieData.getProcessIndex('MaskRefinementProcess');
    if ~isempty(iMaskProcess)
        maskProc = movieData.getProcess(iMaskProcess);
        existMask = true;
    end
end

if ~isempty(iMaskProcess)
    iChan=find(maskProc.checkChannelOutput);
    % if the bead channel is also selected for segmentation, we need to
    % remove that
    iChan = setdiff(iChan,iBeadChan);
    if length(iChan)>1
        iChan=iChan(1);
    end
end

%% Calculate strain energy and total force
tic
minSize = 20; % in pixel
% minTraction = 0; %decided to not use it. used to be 100; % in Pa
forceField=load(forceFieldProc.outFilePaths_{1});
forceField = forceField.forceField;
gridSapcing=forceField(1).pos(2,2)-forceField(1).pos(1,2);
borderWidth=2*gridSapcing;
% nTopBlobs=150; % the number of top maxForceBlobs in single frames
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];

bandwidthNA_pix = round(p.bandWidth*1000/movieData.pixelSize_);

useFOV = p.useFOV;
useCellMask = p.useCellMask;
performForceBlobAnalysis=p.performForceBlobAnalysis;

logMsg='Quantifying strain energy and total force';
if feature('ShowFigureWindows'), waitbar(0,wtBar,sprintf(logMsg)); end
for ii=1:nFrames
    % Make sure if each tmap has its contents
%     curTMap=tMap(:,:,ii);
    curTMap=forceFieldProc.loadChannelOutput('output','tMapUnshifted','iFrame',ii,'noStackRequired',true);

%     curDMap=dMap(:,:,ii);
    curDMap=dispProc.loadChannelOutput('output','dMapUnshifted','iFrame',ii,'noStackRequired',true);
    if isempty(curTMap)
        try
            curTMap=load(outFileTMap(ii),'cur_tMap');
            curTMap = curTMap.cur_tMap;
        catch
            curDisplField = CorrectedDisplFieldProc.loadChannelOutput(ii);
            [curTMapInsert,~,~,cropInfo]=generateHeatmapShifted(forceField(ii),curDisplField,0);
            curTMap = zeros(size(curDMap));
            curTMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3))=curTMapInsert{1};
        end
    end
    if useFOV || 1 % I will calculate this anyway
        % segment
        maskFOV = curTMap>0;%~isnan(curTMap); 
        % To crop, you can calculate the top, bottom, left, and right columns of the mask's "true" area by using any(). 
        anyCol=any(maskFOV);
        anyRow=any(maskFOV,2);
        row1 = find(anyRow,1);
        row2 = find(anyRow,1,'last');
        col1 = find(anyCol,1);
        col2 = find(anyCol,1,'last');

        tMapFOV = curTMap(row1:row2, col1:col2);
        dMapFOV = curDMap(row1:row2, col1:col2);

        % Erode FOV mask to remove the edge effect 
        maskShrunkenBorder = true(size(tMapFOV));
        maskShrunkenBorder = bwmorph(maskShrunkenBorder,'erode',borderWidth);

        curSE_FOV_SE=1/2*sum(sum(dMapFOV.*tMapFOV.*maskShrunkenBorder))*(pixSize_mu*1e-6)^3*1e15; % in femto-Joule=1e15*(N*m)
        curSE_FOV_area=sum(maskShrunkenBorder(:))*areaConvert; % this is in um2
        curSE_FOV_SEDensity=curSE_FOV_SE/curSE_FOV_area*1e3; % J/m2
        
        SE_FOV_SE(ii)=curSE_FOV_SE; % in femto-Joule=1e15*(N*m)
        SE_FOV_area(ii)=curSE_FOV_area; % this is in um2
        SE_FOV_SEDensity(ii)=curSE_FOV_SEDensity; % J/m2
        
        reducedTmap = tMapFOV.*maskShrunkenBorder;
        totalForceFOV(ii) = sum(sum(reducedTmap))*areaConvert*1e-3; % in nN
        avgTractionFOV(ii) = mean(tMapFOV(maskShrunkenBorder)); % in Pa
    else
        SE_FOV_SE(ii)=NaN; % in femto-Joule=1e15*(N*m)
        SE_FOV_area(ii)=NaN; % this is in um2
        SE_FOV_SEDensity(ii)=NaN; % J/m2
        totalForceFOV(ii) = NaN;
        avgTractionFOV(ii) = NaN;
    end
    
    if existMask && useCellMask
        maskCell = maskProc.loadChannelOutput(iChan,ii);

        ref_obj = imref2d(size(maskCell));
        if existSDC
%             if isa(SDCProc,'EfficientSubpixelRegistrationProcess')
%                 maxX = 0;
%                 maxY = 0;
%             else        
%                 maxX = ceil(max(abs(T(:, 2))));
%                 maxY = ceil(max(abs(T(:, 1))));
%             end
            Tr = affine2d([1 0 0; 0 1 0; fliplr(T(ii, :)) 1]);
            % Apply subpixel-wise registration to original masks

%             Ibw = padarray(maskCell, [maxY, maxX]);
            maskCell = imwarp(maskCell, Tr,'OutputView',ref_obj);
        end
        % mask for band from edge
        iMask = imcomplement(maskCell);
        distFromEdge = bwdist(iMask);
        bandMask = distFromEdge <= bandwidthNA_pix;

        maskOnlyBand = bandMask & maskCell;
        maskInterior = ~bandMask & maskCell;
        % saving it
        save(maskIndPath(ii), 'maskOnlyBand')
        
        areaCell = sum(maskCell(:))*areaConvert;  % in um2
        areaPeri = sum(maskOnlyBand(:))*areaConvert; % in um2
        areaInside = sum(maskInterior(:))*areaConvert; % in um2

        tMapCell = curTMap(maskCell);
        dMapCell = curDMap(maskCell);
        tMapCellPeri = curTMap(maskOnlyBand);
        dMapCellPeri = curDMap(maskOnlyBand);
        tMapCellInside = curTMap(maskInterior);
        dMapCellInside = curDMap(maskInterior);

        SE_Cell_SE(ii)=1/2*sum(sum(dMapCell.*tMapCell))*(pixSize_mu*1e-6)^3*1e15; % in femto-Joule=1e15*(N*m)
        SE_Cell_area(ii)=areaCell; % this is in um2
        SE_Cell_SEDensity(ii)=SE_Cell_SE(ii)/SE_Cell_area(ii)*1e3; % J/m2
        SE_Cell_SE_peri(ii)=1/2*sum(sum(dMapCellPeri.*tMapCellPeri))*(pixSize_mu*1e-6)^3*1e15; % in femto-Joule=1e15*(N*m)
        SE_Cell_SE_inside(ii)=1/2*sum(sum(dMapCellInside.*tMapCellInside))*(pixSize_mu*1e-6)^3*1e15; % in femto-Joule=1e15*(N*m)
        SE_Cell_SEDensityPeri(ii)=SE_Cell_SE_peri(ii)/areaPeri*1e3; % J/m2
        SE_Cell_SEDensityInside(ii)=SE_Cell_SE_inside(ii)/areaInside*1e3; % J/m2
        totalForceCell(ii) = sum(sum(tMapCell))*areaConvert*1e-3; % in nN
        totalForceCellPeri(ii) = sum(sum(tMapCellPeri))*areaConvert*1e-3; % in nN
        totalForceCellInside(ii) = sum(sum(tMapCellInside))*areaConvert*1e-3; % in nN

        avgTractionCell(ii) = mean(tMapCell(:)); % in Pa
        avgTractionCellPeri(ii) = mean(tMapCellPeri); % in Pa
        avgTractionCellInside(ii) = mean(tMapCellInside); % in Pa
        
        totalDispCell(ii) = sum(sum(dMapCell))*areaConvert*pixSize_mu; % in um3
        totalDispCellPeri(ii) = sum(sum(dMapCellPeri))*areaConvert*pixSize_mu; % in um3
        totalDispCellInside(ii) = sum(sum(dMapCellInside))*areaConvert*pixSize_mu; % in um3
               
        avgDispCell(ii) = totalDispCell(ii)/areaCell; % um
        avgDispCellPeri(ii) = totalDispCellPeri(ii)/areaPeri; %um
        avgDispCellInside(ii) = totalDispCellInside(ii)/areaInside; %um
    end
        
    if performForceBlobAnalysis
        maskForceBlob = blobSegmentThresholdTFM(tMapFOV,minSize,0,maskShrunkenBorder);
%         maskForceBlob = bwmorph(maskForceBlob,'dilate',1);
%         maskForceBlob = padarray(maskForceBlob,[borderWidth borderWidth]);
%         maskHighTraction=tMapFOV>minTraction;
%         maskForceBlob = maskForceBlob & maskHighTraction;

        SE_Blobs_SE(ii)=1/2*sum(dMapFOV(maskForceBlob).*tMapFOV(maskForceBlob))*(pixSize_mu*1e-6)^3; % this is in Newton*m.
        SE_Blobs_SE(ii)=SE_Blobs_SE(ii)*1e15; % this is now in femto-Joule
        
        stats=regionprops(maskForceBlob,tMapFOV,'Area','PixelIdxList','Centroid','MinIntensity','MaxIntensity','MeanIntensity','WeightedCentroid');
        SE_Blobs_nFA(ii)=numel(stats);
        SE_Blobs_areaFA(ii) = sum(maskForceBlob(:))*areaConvert; % in um2
        SE_Blobs_avgFAarea(ii) = SE_Blobs_areaFA(ii)/SE_Blobs_nFA(ii);
        SE_Blobs_avgSEperFA(ii) = SE_Blobs_SE(ii)/SE_Blobs_nFA(ii); % still in femto-J
        SE_Blobs_SEDensity(ii)=SE_Blobs_SE(ii)/SE_Blobs_areaFA(ii)*1e3; % J/m2
        individualForceBlobs = arrayfun(@(x) x.MeanIntensity,stats);
        individualForceBlobMax = arrayfun(@(x) x.MaxIntensity,stats);
    %     individualForceBlobMin = arrayfun(@(x) x.MinIntensity,stats);
%         individualForceBlobCenters = arrayfun(@(x) x.WeightedCentroid,stats,'UniformOutput',false);
%         individualForceBlobAreas = arrayfun(@(x) x.Area,stats);

        totalForceBlobs_force(ii,1)=sum(tMapFOV(maskForceBlob))*areaConvert*1e-3; % in nN %*(pixSize_mu*1e-6)^2*1e9; % in nN
        totalForceBlobs_avgTraction{ii,1}=(individualForceBlobs); % in Pa
        totalForceBlobs_maxTraction{ii,1}=(individualForceBlobMax); % in Pa
        totalForceBlobs_forceBlobPixelIdxList{ii,1}=arrayfun(@(x) x.PixelIdxList,stats,'UniformOutput',false);
        % find an adhesion that contains top three max traction
    %     [individualForceBlobMaxSorted, topIDs]=sort(individualForceBlobMax,'descend');
        if existMask && useCellMask
            % This is to quantify force blob inside the cell
            maskCellFOV = maskCell(row1:row2, col1:col2);
            maskForceBlobCell = maskForceBlob & maskCellFOV;
            statsCell = regionprops(maskForceBlobCell,tMapFOV,'Area','PixelIdxList','Centroid','MinIntensity','MaxIntensity','MeanIntensity','WeightedCentroid');
            individualForceBlobsCell = arrayfun(@(x) x.MeanIntensity,statsCell);
            totalForceBlobs_avgTractionCell{ii,1}=(individualForceBlobsCell); % in Pa
        end
    end
    
    % See if there is overlap between interiorMask and maskAdhesion - for
    % this I have to use AND operator
    
%     for k=1:SE_Blobs.nFA(ii)
% %         SE_Blobs.individualSE(ii).SE(k)=1/2*sum(dMapFOV(stats(k).PixelIdxList).*tMapFOV(stats(k).PixelIdxList))*(pixSize_mu*1e-6)^3*1e15;
% %         SE_Blobs.individualSE(ii).area(k)=stats(k).Area*areaConvert;
% %         SE_Blobs.individualSE(ii).SED(k)=SE_Blobs.individualSE(ii).SE(k)/SE_Blobs.individualSE(ii).area(k)*1e3;% J/m2.
%     end
%     [maxT]=max(tMapFOV(:));
%     SE_Blobs.maxT(ii) = maxT;
%     [SE_Blobs.maxSEFA(ii),~]=max(SE_Blobs.individualSE(ii).SE);
%     [SE_Blobs.maxSED_FA(ii),~]=max(SE_Blobs.individualSE(ii).SED);
    
    % Update the waitbar
    if feature('ShowFigureWindows')
        tj=toc;
        waitbar(ii/nFrames,wtBar,sprintf([logMsg timeMsg(tj*(nFrames-ii)/ii)]));
    end
end
SE_FOV.SE = SE_FOV_SE;
SE_FOV.area = SE_FOV_area;
SE_FOV.SEDensity=SE_FOV_SEDensity;
SE_Cell.SE=SE_Cell_SE; % in femto-Joule=1e15*(N*m)
SE_Cell.area=SE_Cell_area; % this is in um2
SE_Cell.SEDensity=SE_Cell_SEDensity; % J/m2
SE_Cell.SE_peri=SE_Cell_SE_peri; % in femto-Joule=1e15*(N*m)
SE_Cell.SE_inside=SE_Cell_SE_inside; % in femto-Joule=1e15*(N*m)
SE_Cell.SEDensityPeri=SE_Cell_SEDensityPeri; % J/m2
SE_Cell.SEDensityInside=SE_Cell_SEDensityInside; % J/m2

SE_Blobs.SE=SE_Blobs_SE; % this is in femto-Joule
SE_Blobs.nFA=SE_Blobs_nFA;
SE_Blobs.areaFA=SE_Blobs_areaFA; % in um2
SE_Blobs.avgFAarea=SE_Blobs_avgFAarea;
SE_Blobs.avgSEperFA=SE_Blobs_avgSEperFA; % still in femto-J
SE_Blobs.SEDensity=SE_Blobs_SEDensity; % J/m2

totalForceBlobs.force=        totalForceBlobs_force; % in nN %*(pixSize_mu*1e-6)^2*1e9; % in nN
totalForceBlobs.avgTraction =       totalForceBlobs_avgTraction; % in Pa
totalForceBlobs.maxTraction =         totalForceBlobs_maxTraction; % in Pa
totalForceBlobs.forceBlobPixelIdxList =        totalForceBlobs_forceBlobPixelIdxList;
totalForceBlobs.avgTractionCell =        totalForceBlobs_avgTractionCell; % in Pa

%% Save
logMsg='Saving...';
if feature('ShowFigureWindows'), waitbar(0,wtBar,sprintf(logMsg)); end
if useFOV || 1
    save(outputFile{1},'SE_FOV','totalForceFOV','avgTractionFOV');
    if p.exportCSV
        tableSE_FOV=struct2table(SE_FOV);
        writetable(tableSE_FOV,outputFile{4})
%         writetable(tableForceFOV,outputFile{5})
        totalAvgTractionFOV=cell2table({totalForceFOV, avgTractionFOV},...
            'VariableNames',{'totalForceFOV','avgTractionFOV'});        
        writetable(totalAvgTractionFOV,outputFile{5})
    end
end
if existMask && useCellMask
    save(outputFile{2},'SE_Cell','totalForceCell','totalForceCellPeri',...
        'totalForceCellInside','totalDispCell','totalDispCellPeri','totalDispCellInside',...
        'avgTractionCell','avgTractionCellPeri','avgTractionCellInside',...
        'avgDispCell','avgDispCellPeri','avgDispCellInside'); % need to be updated for faster loading. SH 20141106
    if p.exportCSV
        tableSE_Cell=struct2table(SE_Cell);
        writetable(tableSE_Cell,outputFile{6})
        tableForceCell=table(totalForceCell,'VariableNames',{'totalForceCell'});
        writetable(tableForceCell,outputFile{7})
        tableForceCell2=table(totalForceCellPeri,'VariableNames',{'totalForceCellPeri'});
        writetable(tableForceCell2,[p.OutputDirectory filesep 'totalForceCellPeri.csv'])
        tableForceCell3=table(totalForceCellInside,'VariableNames',{'totalForceCellInside'});
        writetable(tableForceCell3,[p.OutputDirectory filesep 'totalForceCellInside.csv'])
        tableForceCell4=table(avgTractionCell,'VariableNames',{'avgTractionCell'});
        writetable(tableForceCell4,[p.OutputDirectory filesep 'avgTractionCell.csv'])
        tableForceCell5=table(avgTractionCellPeri,'VariableNames',{'avgTractionCellPeri'});
        writetable(tableForceCell5,[p.OutputDirectory filesep 'avgTractionCellPeri.csv'])
        tableForceCell6=table(avgTractionCellInside,'VariableNames',{'avgTractionCellInside'});
        writetable(tableForceCell6,[p.OutputDirectory filesep 'avgTractionCellInside.csv'])
    end
end
if performForceBlobAnalysis
    save(outputFile{3},'SE_Blobs','totalForceBlobs');
    if p.exportCSV
        tableSE_Blobs=struct2table(SE_Blobs);
        writetable(tableSE_Blobs,outputFile{8})
        tableForceBlobs=table(totalForceBlobs.force,'VariableNames',{'totalForceBlobs'});
        writetable(tableForceBlobs,outputFile{9})
    end
end
save(outputFile{1,10},'SE_Blobs','totalForceBlobs', 'SE_Cell','totalForceCell','SE_FOV',...
    'totalForceFOV','totalForceCellPeri','totalForceCellInside','avgTractionCell',...
    'avgTractionCellPeri','avgTractionCellInside','-v7.3')
%% Close waitbar
if feature('ShowFigureWindows'), close(wtBar); end

disp('Finished calculating strain energy and total force!')
