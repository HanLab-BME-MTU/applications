function analysisStruct = makiMetaPlateAnalysis(jobType,analysisStruct,...
    verbose)
%MAKIMETAPLATEANALYSIS analyzes the behavior of the metaphase plate
%
%SYNOPSIS analysisStruct = makiMetaPlateAnalysis(jobType,analysisStruct,...
%    verbose)
%
%INPUT  jobType: string which can take the values:
%               'TEST', 'HERCULES', 'DANUSER', 'MERALDI', 'SWEDLOW' or
%               'MCAINSH'
%       analysisStruct: Structure with fields
%               .fileName: Name of file where analysisStruct is stored.
%               .filePath: Directory where analysisStruct is stored.
%               .movies  : nx2 cell array of the names of the n selected movies.
%                          First column: file name, second column: file path.
%                       Optional. If not input, GUI to load movies is launched.
%       verbose: 1 to make plots, 0 otherwise. Optional. Default: 0.
%
%OUTPUT analysisStruct: Same as input but with additional field
%           .metaPlate: 
%
%REMARKS Code is not applicable to anaphase movies/frames
%
%Khuloud Jaqaman, August 2007

%% input
if nargin < 1 || isempty(jobType)
    jobType = 'DANUSER';
end

if nargin < 3 || isempty(verbose)
    verbose = 0;
end

%interactively obtain analysisStruct if not input
if nargin < 2 || isempty(analysisStruct) || ~isfield(analysisStruct,'movies')
    analysisStruct = makiCollectMovies(jobType);
end

%convert file paths in analysisStruct from identifier to real path
analysisStruct = makiMakeAnalysisPlatformIndependent(analysisStruct,jobType);

%extract fileName and directory name from analysisStruct
fileName = analysisStruct.fileName;
dir2SaveRes = analysisStruct.filePath;

%load dataStruct's belonging to the movies in analysisStruct
moviesList = analysisStruct.movies;
numMovies = size(moviesList,1);
for iMovie = numMovies : -1 : 1
    dataStruct(iMovie,1) = makiLoadDataFile(jobType,fullfile(moviesList{iMovie,2},moviesList{iMovie,1}));
end

%% get metaphase plate information

%reserve memory
metaPlateOrigin = repmat(struct('observations',[]),numMovies,1);
metaPlateNormal = repmat(struct('observations',[]),numMovies,1);
kinetScatterStd = repmat(struct('observations',[]),numMovies,1);

%go over all movies
for iMovie = 1 : numMovies
    
    %copy fields out of dataStruct(iMovie)
    planeFit = dataStruct(iMovie).planeFit;
    updatedClass = dataStruct(iMovie).updatedClass;

    %get number of frames in movie
    numFrames = length(planeFit);
    
    %get frame phases and identify non-anaphase frames that have a plane
    framePhase = vertcat(planeFit.phase);
    framesWithPlane = find(framePhase == 'p' | framePhase == 'm');

    %go over frames with plane and collect information
    planeOrigin = NaN(numFrames,3);
    planeNormal = NaN(numFrames,3);
    kinetStd = NaN(numFrames,1);
    for iFrame = framesWithPlane'
        
        %get plane origin
        planeOrigin(iFrame,:) = planeFit(iFrame).planeOrigin;
        
        %get vector normal to plane
        planeNormal(iFrame,:) = planeFit(iFrame).planeVectors(:,1)';
        
        %get coordinates of inlier kinetochores
        inlierIndx = updatedClass(iFrame).inlierIdx;
        kinetCoord = planeFit(iFrame).planeCoord(inlierIndx,1);
        kinetStd(iFrame) = std(kinetCoord);
        
    end
    
    %store information in global arrays
    metaPlateOrigin(iMovie).observations = planeOrigin;
    metaPlateNormal(iMovie).observations = planeNormal;
    kinetScatterStd(iMovie).observations = kinetStd;
    
end
    
%% origin displacement and normal orientation change

%reserve memory
originDispMag = repmat(struct('observations',[]),numMovies,1);
originDispAlongNorm = repmat(struct('observations',[]),numMovies,1);
originDispPerpNorm = repmat(struct('observations',[]),numMovies,1);
normalDirChange = repmat(struct('observations',[]),numMovies,1);

%go over all movies
for iMovie = 1 : numMovies
   
    %get the frame-to-frame displacement of the plane origin
    originDisp = metaPlateOrigin(iMovie).observations(2:end,:) - ...
        metaPlateOrigin(iMovie).observations(1:end-1,:);
    
    %get the normal direction at every frame
    normalDir = metaPlateNormal(iMovie).observations;
    
    %get the displacement magnitude
    originDispMagT = sqrt(sum(originDisp.^2,2)); %um
    
    %get the component of the displacement along the normal
    originDispAlongNormT = sum(originDisp .* normalDir(1:end-1,:),2); %um
    
    %get the component of the displacement perpendicular to the normal
    originDispPerpNormT = sqrt(originDispMagT.^2 - originDispAlongNormT.^2); %um
    
    %calculate the angle between the normals in consecutive frames
    normalDirChangeT = acos( sum( normalDir(1:end-1,:) .* ...
        normalDir(2:end,:), 2) ) * 180 / pi; %degrees
    
    %store data
    originDispMag(iMovie).observations = originDispMagT;
    originDispAlongNorm(iMovie).observations = originDispAlongNormT;
    originDispPerpNorm(iMovie).observations = originDispPerpNormT;
    normalDirChange(iMovie).observations = normalDirChangeT;
        
end

%% distributions and some distribution parameters

%initialization
originDispMagIndParam = NaN(numMovies,7);
originDispAlongNormIndParam = NaN(numMovies,7);
originDispPerpNormIndParam = NaN(numMovies,7);
normalDirChangeIndParam = NaN(numMovies,7);
kinetScatterStdIndParam = NaN(numMovies,7);

% overall %

%displacement magnitude
originDispMagDistr = vertcat(originDispMag.observations);
originDispMagDistr = originDispMagDistr(~isnan(originDispMagDistr));
originDispMagParam = [mean(originDispMagDistr) std(originDispMagDistr) ...
    min(originDispMagDistr) prctile(originDispMagDistr,25) ...
    prctile(originDispMagDistr,50) prctile(originDispMagDistr,75) ...
    max(originDispMagDistr)];

%displacement along normal
originDispAlongNormDistr = vertcat(originDispAlongNorm.observations);
originDispAlongNormDistr = originDispAlongNormDistr(~isnan(originDispAlongNormDistr));
originDispAlongNormParam = [mean(originDispAlongNormDistr) std(originDispAlongNormDistr) ...
    min(originDispAlongNormDistr) prctile(originDispAlongNormDistr,25) ...
    prctile(originDispAlongNormDistr,50) prctile(originDispAlongNormDistr,75) ...
    max(originDispAlongNormDistr)];

%displacement perpendicular to normal
originDispPerpNormDistr = vertcat(originDispPerpNorm.observations);
originDispPerpNormDistr = originDispPerpNormDistr(~isnan(originDispPerpNormDistr));
originDispPerpNormParam = [mean(originDispPerpNormDistr) std(originDispPerpNormDistr) ...
    min(originDispPerpNormDistr) prctile(originDispPerpNormDistr,25) ...
    prctile(originDispPerpNormDistr,50) prctile(originDispPerpNormDistr,75) ...
    max(originDispPerpNormDistr)];

%normal direction change
normalDirChangeDistr = vertcat(normalDirChange.observations);
normalDirChangeDistr = normalDirChangeDistr(~isnan(normalDirChangeDistr));
normalDirChangeParam = [mean(normalDirChangeDistr) std(normalDirChangeDistr) ...
    min(normalDirChangeDistr) prctile(normalDirChangeDistr,25) ...
    prctile(normalDirChangeDistr,50) prctile(normalDirChangeDistr,75) ...
    max(normalDirChangeDistr)];

%kinetochore scatter std
kinetScatterStdDistr = vertcat(kinetScatterStd.observations);
kinetScatterStdDistr = kinetScatterStdDistr(~isnan(kinetScatterStdDistr));
kinetScatterStdParam = [mean(kinetScatterStdDistr) std(kinetScatterStdDistr) ...
    min(kinetScatterStdDistr) prctile(kinetScatterStdDistr,25) ...
    prctile(kinetScatterStdDistr,50) prctile(kinetScatterStdDistr,75) ...
    max(kinetScatterStdDistr)];

% individual cells %

%displacement magnitude
for iMovie = 1 : numMovies
    allValues = originDispMag(iMovie).observations;
    allValues = allValues(~isnan(allValues));
    if ~isempty(allValues)
        originDispMagIndParam(iMovie,:) = [mean(allValues) std(allValues) ...
            min(allValues) prctile(allValues,25) prctile(allValues,50) ...
            prctile(allValues,75) max(allValues)];
    end
end

%displacement along normal
for iMovie = 1 : numMovies
    allValues = originDispAlongNorm(iMovie).observations;
    allValues = allValues(~isnan(allValues));
    if ~isempty(allValues)
        originDispAlongNormIndParam(iMovie,:) = [mean(allValues) std(allValues) ...
            min(allValues) prctile(allValues,25) prctile(allValues,50) ...
            prctile(allValues,75) max(allValues)];
    end
end

%displacement perpendicular to normal
for iMovie = 1 : numMovies
    allValues = originDispPerpNorm(iMovie).observations;
    allValues = allValues(~isnan(allValues));
    if ~isempty(allValues)
        originDispPerpNormIndParam(iMovie,:) = [mean(allValues) std(allValues) ...
            min(allValues) prctile(allValues,25) prctile(allValues,50) ...
            prctile(allValues,75) max(allValues)];
    end
end

%normal direction change
for iMovie = 1 : numMovies
    allValues = normalDirChange(iMovie).observations;
    allValues = allValues(~isnan(allValues));
    if ~isempty(allValues)
        normalDirChangeIndParam(iMovie,:) = [mean(allValues) std(allValues) ...
            min(allValues) prctile(allValues,25) prctile(allValues,50) ...
            prctile(allValues,75) max(allValues)];
    end
end

%kinetochore scatter std
for iMovie = 1 : numMovies
    allValues = kinetScatterStd(iMovie).observations;
    allValues = allValues(~isnan(allValues));
    if ~isempty(allValues)
        kinetScatterStdIndParam(iMovie,:) = [mean(allValues) std(allValues) ...
            min(allValues) prctile(allValues,25) prctile(allValues,50) ...
            prctile(allValues,75) max(allValues)];
    end
end

%% autocorrelation of origin displacement along normal

%define maximum lag
maxLag = 20;

%displacement along normal autocorrelation
[originDispAlongNormAutocorr,errFlag] = autoCorr(originDispAlongNorm,maxLag);

%% output to analysisStruct

distribution = struct('originDispMag',originDispMagDistr,...
    'originDispAlongNorm',originDispAlongNormDistr,...
    'originDispPerpNorm',originDispPerpNormDistr,...
    'normalDirChange',normalDirChangeDistr,...
    'kinetScatterStd',kinetScatterStdDistr);
meanStdMin25P50P75PMax.all = struct('originDispMag',originDispMagParam,...
    'originDispAlongNorm',originDispAlongNormParam,...
    'originDispPerpNorm',originDispPerpNormParam,...
    'normalDirChange',normalDirChangeParam,...
    'kinetScatterStd',kinetScatterStdParam);
meanStdMin25P50P75PMax.indcell = struct('originDispMag',originDispMagIndParam,...
    'originDispAlongNorm',originDispAlongNormIndParam,...
    'originDispPerpNorm',originDispPerpNormIndParam,...
    'normalDirChange',normalDirChangeIndParam,...
    'kinetScatterStd',kinetScatterStdIndParam);
autocorr = struct('originDispAlongNorm',originDispAlongNormAutocorr);

metaPlate = struct('distribution',distribution,...
    'meanStdMin25P50P75PMax',meanStdMin25P50P75PMax,...
    'autocorr',autocorr);

%check whether current analysisStruct already has the matePlate field
fieldExists = isfield(analysisStruct,'metaPlate');

%store results in analysisStruct
analysisStruct.metaPlate = metaPlate;

%if metaPlate field already existed, add 1 to the version number in
%the file name where analysisStruct will be stored
if fieldExists
    [versionNum,fileBody] = makiGetVersion(fileName);
    fileName = [fileBody '_' num2str(versionNum+1) '.mat'];
    analysisStruct.fileName = fileName;
end

%save analysisStruct
analysisStruct = makiMakeAnalysisPlatformIndependent(analysisStruct,jobType);
save(fullfile(dir2SaveRes,fileName),'analysisStruct');

%% plots

if verbose

    %get time between frames
    timeLapse = round(2*dataStruct(1).dataProperties.timeLapse)/2;

    %get number of frames in each movie
    numFrames = NaN(numMovies,1);
    for iMovie = 1 : numMovies
        numFrames(iMovie) = dataStruct(iMovie).dataProperties.movieSize(end);
    end
    numFrames = min(numFrames)-1;

    %open figure and write title
    figFileName = [fileName(1:end-4) '-PlateDisplacement'];
    figHandle = figure('Name',figFileName,'NumberTitle','off');

        %create subplot 1
        subplot(2,2,1);
        hold on

        %put all displacement magnitudes together in one matrix
        dispMat = [];
        for iMovie = 1 : numMovies
            dispMat = [dispMat originDispMag(iMovie).observations(1:numFrames,1)];
        end

        %plot distance over time for all movies
        plot((0:numFrames-1)*timeLapse,dispMat);

        %set axes limit
        axis([0 (numFrames-1)*timeLapse 0 max(dispMat(:))+1]);

        %write axes labels
        xlabel('Time (s)');
        ylabel('Plate displacement magnitude (\mum)');

        %hold off subplot 1
        hold off

        %create subplot 2
        subplot(2,2,2);
        hold on

        %put all displacements along normal together in one matrix
        dispMat = [];
        for iMovie = 1 : numMovies
            dispMat = [dispMat originDispAlongNorm(iMovie).observations(1:numFrames,1)];
        end

        %plot distance over time for all movies
        plot((0:numFrames-1)*timeLapse,dispMat);

        %set axes limit
        axis([0 (numFrames-1)*timeLapse min(dispMat(:))-1 max(dispMat(:))+1]);

        %write axes labels
        xlabel('Time (s)');
        ylabel('Plate displacement along normal (\mum)');

        %hold off subplot 2
        hold off

        %create subplot 3
        subplot(2,2,3);
        hold on

        %put all normal direction changes together in one matrix
        dispMat = [];
        for iMovie = 1 : numMovies
            dispMat = [dispMat normalDirChange(iMovie).observations(1:numFrames,1)];
        end

        %plot distance over time for all movies
        plot((0:numFrames-1)*timeLapse,dispMat);

        %set axes limit
        axis([0 (numFrames-1)*timeLapse 0 max(dispMat(:))+1]);

        %write axes labels
        xlabel('Time (s)');
        ylabel('Plate normal direction change (degrees)');

        %hold off subplot 3
        hold off

        %create subplot 4
        subplot(2,2,4);
        hold on

        %plot displacement along normal autocorrelation
        if ~isempty(originDispAlongNormAutocorr)
            plot((0:maxLag)*timeLapse,originDispAlongNormAutocorr(:,1));

            %set axes limit
            axis([0 (maxLag)*timeLapse min(0,1.1*min(originDispAlongNormAutocorr(:,1))) 1.1]);

            %write axes labels
            xlabel('Time (s)');
            ylabel('Autocorrelation of plate displacement along normal');
        end
        
        %hold off subplot 4
        hold off

        %save figure in file
        saveas(figHandle,fullfile(dir2SaveRes,figFileName),'fig');

end

%% ~~~ the end ~~~ %%

