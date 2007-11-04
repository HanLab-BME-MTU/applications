function analysisStruct = makiDiffusionAnalysis(jobType,analysisStruct,verbose)
%MAKIDIFFUSIONANALYSIS calculate the diffusion parameters for kinetochore motion
%
%SYNOPSIS analysisStruct = makiDiffusionAnalysis(jobType,analysisStruct,verbose)
%
%INPUT  jobType: string which can take the values:
%               'TEST','HERCULES','DANUSER','MERLADI',SWEDLOW' or MCAINSH
%       analysisStruct: Structure with fields
%               .fileName: Name of file where analysisStruct is stored.
%               .filePath: Directory where analysisStruct is stored.
%               .movies  : nx2 cell array of the names of the n selected movies.
%                          First column: file name, second column: file path.
%                       Optional. If not input, GUI to load movies is launched.
%       verbose: 1 to make plots, 0 otherwise. Optional. Default: 0.
%
%OUTPUT analysisStruct: Same as input but with additional field
%           .diffusionAnalysis: 
%
%Khuloud Jaqaman, November 2007

%% input
if nargin < 1 || isempty(jobType)
    jobType = 'TEST';
end

if nargin < 3 || isempty(verbose)
    verbose = 0;
end

%interactively obtain analysisStruct if not input
if nargin < 2 || isempty(analysisStruct) || ~isfield(analysisStruct,'movies')
    analysisStruct = makiCollectMovies(jobType);
end

%convert file paths in analysisStruct from identifier to real path
analysisStruct = makiMakeAnalysisPlatformIndependent(analysisStruct);

%extract fileName and directory name from analysisStruct
fileName = analysisStruct.fileName;
dir2SaveRes = analysisStruct.filePath;

%load dataStruct's belonging to the movies in analysisStruct
moviesList = analysisStruct.movies;
numMovies = size(moviesList,1);
for iMovie = numMovies : -1 : 1
    dataStruct(iMovie,1) = makiLoadDataFile(fullfile(moviesList{iMovie,2},moviesList{iMovie,1}));
end

%find number of tracks in each movie
numTracks = zeros(numMovies,1);
for iMovie = 1 : numMovies
    numTracks(iMovie) = length(dataStruct(iMovie).tracks);
end
numTracksTot = sum(numTracks);

%find number of frames in each movie
numFrames = zeros(numMovies,1);
for iMovie = 1 : numMovies
    numFrames = dataStruct(iMovie).dataProperties.movieSize(end);
end
numFramesMax = max(numFrames);

%get time between frames
timeLapse = round(dataStruct(1).dataProperties.timeLapse);

%% collect tracks

%initialize global track index
iGlobalInlier = 0;
iGlobalUnaligned = 0;
iGlobalLagging = 0;

%reserve memory for big track matrices
tracksInlier = NaN(numTracksTot,numFramesMax*8);
tracksUnaligned = NaN(numTracksTot,numFramesMax*8);
tracksLagging = NaN(numTracksTot,numFramesMax*8);

%go over all movies
for iMovie = 1 : numMovies

    %copy fields out of dataStruct(iMovie)
    tracks = dataStruct(iMovie).tracks;
    frameAlignment = dataStruct(iMovie).frameAlignment;
    updatedClass = dataStruct(iMovie).updatedClass;

    if numTracks(iMovie) > 0

        %track type = 0 if all kinetochores are inliers
        %track type = 1 if some kinetochores are unaligned
        %track type = 2 if some kinetochores are lagging
        trackType = zeros(numTracks(iMovie),1);
        trackType(updatedClass(1).tracksUnaligned) = 1;
        trackType(updatedClass(1).tracksLagging) = 2;
        
        %get track start, end and life times
        trackSEL = getTrackSEL(tracks);
        
        %convert tracks into matrix format
        [tracksMat,tracksIndxMat] = convStruct2MatNoMS(tracks);
        
        %go over all tracks in movie
        for iTrack = 1 : numTracks(iMovie)

            %find frames where track "exists"
            goodFrames = find(tracksIndxMat(iTrack,:) ~= 0);

            %get aligned coordinates
            coords = NaN(1,trackSEL(iTrack,3)*8);
            for iFrame = 1 : length(goodFrames)
                jFrame = goodFrames(iFrame);
                coords((iFrame-1)*8+1:(iFrame-1)*8+3) = frameAlignment(jFrame).alignedCoord(tracksIndxMat(iTrack,jFrame),1:3);
                coords((iFrame-1)*8+5:(iFrame-1)*8+7) = frameAlignment(jFrame).alignedCoord(tracksIndxMat(iTrack,jFrame),4:6);
            end

            %store tracks based on track type
            switch trackType(iTrack)

                case 0

                    %increase global index of inlier sisters by 1
                    iGlobalInlier = iGlobalInlier + 1;

                    %store information
                    tracksInlier(iGlobalInlier,1:trackSEL(iTrack,3)*8) = coords; %um

                case 1

                    %increase global index of unaligned sisters by 1
                    iGlobalUnaligned = iGlobalUnaligned + 1;

                    %store information
                    tracksUnaligned(iGlobalUnaligned,1:trackSEL(iTrack,3)*8) = coords; %um

                case 2

                    %increase global index of unaligned sisters by 1
                    iGlobalLagging = iGlobalLagging + 1;

                    %store information
                    tracksLagging(iGlobalLagging,1:trackSEL(iTrack,3)*8) = coords; %um
                    
            end %(switch trackType(iTrack))

        end %(for iTrack = 1 : numTracks(iMovie) )

    end %(if numTracks(iMovie) > 0)

end %(for iMovie = 1 : numMovies)

%define cell array of track labels
label{1,1} = 'Inlier';
label{2,1} = 'Unaligned';
label{3,1} = 'Lagging';
for i = 1 : 3
    eval(['label{i,2} = iGlobal' label{i,1} ' ~= 0;']);
end
goodLabel = find(vertcat(label{:,2}))';

%remove unused rows from track matrices
for iLabel = 1 : 3
    eval(['tracks' label{iLabel,1} ' = tracks' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ',:);'])
end
    
%% calculate mean square displacement over time

%define maximum lag
maxLag = numFramesMax-1;

%initialization
for iLabel = 1 : 3
    eval(['meanSqDispEnsemble' label{iLabel,1} ' = [];'])
    eval(['meanSqDispTime' label{iLabel,1} ' = [];'])
    eval(['diffCoefEnsemble' label{iLabel,1} ' = [];'])
    eval(['diffCoefTime' label{iLabel,1} ' = [];'])
end

%calculation
for iLabel = goodLabel

    %ensemble averaging
    eval(['meanSqDispEnsemble' label{iLabel,1} ' = getAllTracksMSqD(tracks' ...
        label{iLabel,1} ',maxLag,1);'])

    %tiem averaging
    eval(['meanSqDispTime' label{iLabel,1} ' = getAllTracksMSqD(tracks' ...
        label{iLabel,1} ',maxLag);'])
    
    %fit first 6 time points with a straight line to estimate diffusion
    %coefficient
    
    %ensemble average
    eval(['[lineParam,S] = polyfit((1:6)''*timeLapse,'...
        'meanSqDispEnsemble' label{iLabel,1} '(1:6,1),1);']);
    varCovMat = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;
    eval(['diffCoefEnsemble' label{iLabel,1} ...
        '.value = [lineParam(1) sqrt(varCovMat(1))];'])
    eval(['diffCoefEnsemble' label{iLabel,1} ...
        '.fitResults = struct(''lineParam'',lineParam,''S'',S);'])
    
    %time average
    eval(['[lineParam,S] = polyfit((1:6)''*timeLapse,'...
        'meanSqDispTime' label{iLabel,1} '(1:6,1),1);']);
    varCovMat = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;
    eval(['diffCoefTime' label{iLabel,1} ...
        '.value = [lineParam(1) sqrt(varCovMat(1))];'])
    eval(['diffCoefTime' label{iLabel,1} ...
        '.fitResults = struct(''lineParam'',lineParam,''S'',S);'])
    
end

%% extract diffusion parameters


%% output to analysisStruct

for iLabel = 1 : 3

    %store results in one structure
    eval(['ensembleAve = struct(''meanSqDisp'',meanSqDispEnsemble' label{iLabel,1} ','...
        '''diffusionParam'',diffCoefEnsemble' label{iLabel,1} ');']);
    eval(['timeAve = struct(''meanSqDisp'',meanSqDispTime' label{iLabel,1} ','...
        '''diffusionParam'',diffCoefTime' label{iLabel,1} ');']);
    eval(['numTracksInd = iGlobal' label{iLabel,1} ';']);

    eval([label{iLabel,1} ' = struct(''numTracks'',numTracksInd,' ...
        '''ensembleAve'',ensembleAve,''timeAve'',timeAve);'])

end

diffusionAnalysis = struct('Inlier',Inlier,'Unaligned',Unaligned,...
    'Lagging',Lagging);

%check whether current analysisStruct already has the sisterConnection field
fieldExists = isfield(analysisStruct,'diffusionAnalysis');

%store results in analysisStruct
analysisStruct.diffusionAnalysis = diffusionAnalysis;

%if diffusionAnalysis field already existed, add 1 to the version number in
%the file name where analysisStruct will be stored
if fieldExists
    [versionNum,fileBody] = makiGetVersion(fileName);
    fileName = [fileBody '_' num2str(versionNum+1) '.mat'];
    analysisStruct.fileName = fileName;
end

%save analysisStruct
analysisStruct = makiMakeAnalysisPlatformIndependent(analysisStruct);
save(fullfile(dir2SaveRes,fileName),'analysisStruct');

%% plots

if verbose

    for iLabel = goodLabel

        %open figure and write title
        figure('Name',[fileName(1:end-4) ' - Mean square displacement - ' label{iLabel,1}],'NumberTitle','off');

        %create subplot 1
        subplot(2,1,1);
        hold on;

        %plot ensemble average
        eval(['ensembleAverage = meanSqDispEnsemble' label{iLabel,1} ';']);
        plot((1:maxLag)*timeLapse,ensembleAverage(:,1),'k');
        myErrorbar((1:maxLag)*timeLapse,ensembleAverage(:,1),...
            ensembleAverage(:,2)./sqrt(ensembleAverage(:,3)));
        
        %plot straight line of pure, unconstrained Brownian motion
        eval(['straightLine = polyval(diffCoefEnsemble' label{iLabel,1} ...
            '.fitResults.lineParam,(1:round(maxLag/2))''*timeLapse);'])
        plot((1:round(maxLag/2))*timeLapse,straightLine,'r--')

        %set axes limit
        axis([0 maxLag*timeLapse 0 1.2*max(ensembleAverage(:,1))]);

        %write axes labels
        xlabel('Lag (s)');
        ylabel('Ensemble average (\mum^2)');

        %write legend
        text(maxLag*timeLapse/2,0.4*max(ensembleAverage(:,1)),...
            sprintf(' Black: Kinetochore motion \n Red: Pure Brownian motion'));

        %hold off figure
        hold off

        %create subplot 2
        subplot(2,1,2);
        hold on;

        %plot time average
        eval(['timeAverage = meanSqDispTime' label{iLabel,1} ';']);
        plot((1:maxLag)*timeLapse,timeAverage(:,1),'k');
        myErrorbar((1:maxLag)*timeLapse,timeAverage(:,1),...
            timeAverage(:,2)./sqrt(timeAverage(:,3)));

        %plot straight line of pure, unconstrained Brownian motion
        eval(['straightLine = polyval(diffCoefTime' label{iLabel,1} ...
            '.fitResults.lineParam,(1:round(maxLag/2))''*timeLapse);'])
        plot((1:round(maxLag/2))*timeLapse,straightLine,'r--')

        %set axes limit
        axis([0 maxLag*timeLapse 0 1.2*max(timeAverage(:,1))]);
        
        %write axes labels
        xlabel('Lag (s)');
        ylabel('Time average (\mum^2)');

        %write legend
        text(maxLag*timeLapse/2,0.4*max(timeAverage(:,1)),...
            sprintf(' Black: Kinetochore motion \n Red: Pure Brownian motion'));

        %hold off figure
        hold off

    end
    
end

%% ~~~ the end ~~~ %%

