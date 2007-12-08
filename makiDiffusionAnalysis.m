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
%REMARKS Code is not applicable to anaphase movies/frames
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

%define cell array of sister labels
label{1,1} = 'Inlier';
label{2,1} = 'Unaligned';
label{3,1} = 'Lagging';

%initialize global sister index and track matrices
for iLabel = 1 : 3

    eval(['iGlobal' label{iLabel,1} ' = 0;'])

    eval(['tracksAlign' label{iLabel,1} ' = NaN(numTracksTot,numFramesMax*8);'])
    eval(['tracksOrig' label{iLabel,1} ' = NaN(numTracksTot,numFramesMax*8);'])
    
end
    
%go over all movies
for iMovie = 1 : numMovies

    %copy fields out of dataStruct(iMovie)
    tracks = dataStruct(iMovie).tracks;
    frameAlignment = dataStruct(iMovie).frameAlignment;
    updatedClass = dataStruct(iMovie).updatedClass;

    %find frame where anaphase starts (if it starts at all)
    framePhase = vertcat(updatedClass.phase);
    firstFrameAna = find(framePhase=='a',1,'first');

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
            coords((firstFrameAna-1)*8+1:end) = NaN;

            %store tracks based on track type
            iLabel = trackType(iTrack) + 1;

            %increase global index of track type by 1
            eval(['iGlobal' label{iLabel,1} ' = iGlobal' label{iLabel,1} ' + 1;'])

            %store information
            eval(['tracksAlign' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                ',1:trackSEL(iTrack,3)*8) = coords;']) %um
            eval(['tracksOrig' label{iLabel,1} '(iGlobal' label{iLabel,1} ...
                ',1:trackSEL(iTrack,3)*8) = tracksMat(iTrack,'...
                '(trackSEL(iTrack,1)-1)*8+1:trackSEL(iTrack,2)*8);'])

        end %(for iTrack = 1 : numTracks(iMovie) )

    end %(if numTracks(iMovie) > 0)

end %(for iMovie = 1 : numMovies)

%store number of tracks per category
for i = 1 : 3
    eval(['label{i,2} = iGlobal' label{i,1} ' ~= 0;']);
end
goodLabel = find(vertcat(label{:,2}))';

%remove unused rows from track matrices
for iLabel = 1 : 3
    eval(['tracksAlign' label{iLabel,1} ' = tracksAlign' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ',:);'])
    eval(['tracksOrig' label{iLabel,1} ' = tracksOrig' label{iLabel,1} ...
        '(1:iGlobal' label{iLabel,1} ',:);'])
end

    
%% mean square displacement & diffusion coefficient

%initialization
for iLabel = 1 : 3
    eval(['meanSqDispAlignEnsemble' label{iLabel,1} ' = [];'])
    eval(['meanSqDispAlignTime' label{iLabel,1} ' = [];'])
    eval(['diffCoefAlignEnsemble' label{iLabel,1} ' = [];'])
    eval(['diffCoefAlignTime' label{iLabel,1} ' = [];'])
    eval(['meanSqDispOrigEnsemble' label{iLabel,1} ' = [];'])
    eval(['meanSqDispOrigTime' label{iLabel,1} ' = [];'])
    eval(['diffCoefOrigEnsemble' label{iLabel,1} ' = [];'])
    eval(['diffCoefOrigTime' label{iLabel,1} ' = [];'])
end

%define maximum lag
maxLag = numFramesMax-1;

%calculation
for iLabel = goodLabel

    %ensemble averaging
    eval(['meanSqDispAlignEnsemble' label{iLabel,1} ' = getAllTracksMSqD(tracksAlign' ...
        label{iLabel,1} ',maxLag,1);'])

    %time averaging
    eval(['meanSqDispAlignTime' label{iLabel,1} ' = getAllTracksMSqD(tracksAlign' ...
        label{iLabel,1} ',maxLag);'])
    
    %ensemble averaging
    eval(['meanSqDispOrigEnsemble' label{iLabel,1} ' = getAllTracksMSqD(tracksOrig' ...
        label{iLabel,1} ',maxLag,1);'])

    %time averaging
    eval(['meanSqDispOrigTime' label{iLabel,1} ' = getAllTracksMSqD(tracksOrig' ...
        label{iLabel,1} ',maxLag);'])
    
    %fit first 6 time points with a straight line to estimate diffusion
    %coefficient
    
    %ensemble average
    eval(['[lineParam,S] = polyfit((1:6)''*timeLapse,'...
        'meanSqDispAlignEnsemble' label{iLabel,1} '(1:6,1),1);']);
    varCovMat = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;
    eval(['diffCoefAlignEnsemble' label{iLabel,1} ...
        '.value = [lineParam(1) sqrt(varCovMat(1))];'])
    eval(['diffCoefAlignEnsemble' label{iLabel,1} ...
        '.fitResults = struct(''lineParam'',lineParam,''S'',S);'])
    eval(['[lineParam,S] = polyfit((1:6)''*timeLapse,'...
        'meanSqDispOrigEnsemble' label{iLabel,1} '(1:6,1),1);']);
    varCovMat = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;
    eval(['diffCoefOrigEnsemble' label{iLabel,1} ...
        '.value = [lineParam(1) sqrt(varCovMat(1))];'])
    eval(['diffCoefOrigEnsemble' label{iLabel,1} ...
        '.fitResults = struct(''lineParam'',lineParam,''S'',S);'])
    
    %time average
    eval(['[lineParam,S] = polyfit((1:6)''*timeLapse,'...
        'meanSqDispAlignTime' label{iLabel,1} '(1:6,1),1);']);
    varCovMat = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;
    eval(['diffCoefAlignTime' label{iLabel,1} ...
        '.value = [lineParam(1) sqrt(varCovMat(1))];'])
    eval(['diffCoefAlignTime' label{iLabel,1} ...
        '.fitResults = struct(''lineParam'',lineParam,''S'',S);'])
    eval(['[lineParam,S] = polyfit((1:6)''*timeLapse,'...
        'meanSqDispOrigTime' label{iLabel,1} '(1:6,1),1);']);
    varCovMat = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;
    eval(['diffCoefOrigTime' label{iLabel,1} ...
        '.value = [lineParam(1) sqrt(varCovMat(1))];'])
    eval(['diffCoefOrigTime' label{iLabel,1} ...
        '.fitResults = struct(''lineParam'',lineParam,''S'',S);'])
    
end

%% output to analysisStruct

for iLabel = 1 : 3

    %store results in one structure
    eval(['numTracksInd = iGlobal' label{iLabel,1} ';']);
    eval(['ensembleAve = struct(''meanSqDisp'',meanSqDispAlignEnsemble' label{iLabel,1} ','...
        '''diffusionParam'',diffCoefAlignEnsemble' label{iLabel,1} ');']);
    eval(['timeAve = struct(''meanSqDisp'',meanSqDispAlignTime' label{iLabel,1} ','...
        '''diffusionParam'',diffCoefAlignTime' label{iLabel,1} ');']);
    alignedCoord = struct('ensembleAve',ensembleAve,'timeAve',timeAve);
    eval(['ensembleAve = struct(''meanSqDisp'',meanSqDispOrigEnsemble' label{iLabel,1} ','...
        '''diffusionParam'',diffCoefOrigEnsemble' label{iLabel,1} ');']);
    eval(['timeAve = struct(''meanSqDisp'',meanSqDispOrigTime' label{iLabel,1} ','...
        '''diffusionParam'',diffCoefOrigTime' label{iLabel,1} ');']);
    originalCoord = struct('ensembleAve',ensembleAve,'timeAve',timeAve);

    eval([label{iLabel,1} ' = struct('...
        '''numTracks'',numTracksInd,' ...
        '''alignedCoord'',alignedCoord,'...
        '''originalCoord'',originalCoord);'])

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
analysisStruct = makiMakeAnalysisPlatformIndependent(analysisStruct,jobType);
save(fullfile(dir2SaveRes,fileName),'analysisStruct');

%% plots

if verbose

    for iLabel = goodLabel

        %open figure and write title
        figure('Name',[fileName(1:end-4) ' - Mean square displacement - ' label{iLabel,1}],'NumberTitle','off');

        %create subplot 1
        subplot(2,2,1);
        hold on;

        %plot ensemble average
        eval(['ensembleAverage = meanSqDispAlignEnsemble' label{iLabel,1} ';']);
        plot((1:maxLag)*timeLapse,ensembleAverage(:,1),'k');
        myErrorbar((1:maxLag)*timeLapse,ensembleAverage(:,1),...
            ensembleAverage(:,2)./sqrt(ensembleAverage(:,3)));
        
        %plot straight line of pure, unconstrained Brownian motion
        eval(['straightLine = polyval(diffCoefAlignEnsemble' label{iLabel,1} ...
            '.fitResults.lineParam,(1:round(maxLag/2))''*timeLapse);'])
        plot((1:round(maxLag/2))*timeLapse,straightLine,'r--')

        %set axes limit
        axis([0 maxLag*timeLapse 0 1.2*max(ensembleAverage(:,1))]);

        %write axes labels
        xlabel('Lag (s)');
        ylabel('Ensemble average (\mum^2)');
        title('From aligned coordinates')

        %write legend
        text(maxLag*timeLapse/2,0.4*max(ensembleAverage(:,1)),...
            sprintf(' Black: Kinetochore motion \n Red: Pure Brownian motion'));

        %hold off figure
        hold off

        %create subplot 2
        subplot(2,2,2);
        hold on;

        %plot ensemble average
        eval(['ensembleAverage = meanSqDispOrigEnsemble' label{iLabel,1} ';']);
        plot((1:maxLag)*timeLapse,ensembleAverage(:,1),'k');
        myErrorbar((1:maxLag)*timeLapse,ensembleAverage(:,1),...
            ensembleAverage(:,2)./sqrt(ensembleAverage(:,3)));
        
        %plot straight line of pure, unconstrained Brownian motion
        eval(['straightLine = polyval(diffCoefOrigEnsemble' label{iLabel,1} ...
            '.fitResults.lineParam,(1:round(maxLag/2))''*timeLapse);'])
        plot((1:round(maxLag/2))*timeLapse,straightLine,'r--')

        %set axes limit
        axis([0 maxLag*timeLapse 0 1.2*max(ensembleAverage(:,1))]);

        %write axes labels
        xlabel('Lag (s)');
        ylabel('Ensemble average (\mum^2)');
        title('From original coordinates')

        %write legend
        text(maxLag*timeLapse/2,0.4*max(ensembleAverage(:,1)),...
            sprintf(' Black: Kinetochore motion \n Red: Pure Brownian motion'));

        %hold off figure
        hold off

        %create subplot 3
        subplot(2,2,3);
        hold on;

        %plot time average
        eval(['timeAverage = meanSqDispAlignTime' label{iLabel,1} ';']);
        plot((1:maxLag)*timeLapse,timeAverage(:,1),'k');
        myErrorbar((1:maxLag)*timeLapse,timeAverage(:,1),...
            timeAverage(:,2)./sqrt(timeAverage(:,3)));

        %plot straight line of pure, unconstrained Brownian motion
        eval(['straightLine = polyval(diffCoefAlignTime' label{iLabel,1} ...
            '.fitResults.lineParam,(1:round(maxLag/2))''*timeLapse);'])
        plot((1:round(maxLag/2))*timeLapse,straightLine,'r--')

        %set axes limit
        axis([0 maxLag*timeLapse 0 1.2*max(timeAverage(:,1))]);
        
        %write axes labels
        xlabel('Lag (s)');
        ylabel('Time average (\mum^2)');
        title('From aligned coordinates')

        %write legend
        text(maxLag*timeLapse/2,0.4*max(timeAverage(:,1)),...
            sprintf(' Black: Kinetochore motion \n Red: Pure Brownian motion'));

        %hold off figure
        hold off

        %create subplot 4
        subplot(2,2,4);
        hold on;

        %plot time average
        eval(['timeAverage = meanSqDispOrigTime' label{iLabel,1} ';']);
        plot((1:maxLag)*timeLapse,timeAverage(:,1),'k');
        myErrorbar((1:maxLag)*timeLapse,timeAverage(:,1),...
            timeAverage(:,2)./sqrt(timeAverage(:,3)));

        %plot straight line of pure, unconstrained Brownian motion
        eval(['straightLine = polyval(diffCoefOrigTime' label{iLabel,1} ...
            '.fitResults.lineParam,(1:round(maxLag/2))''*timeLapse);'])
        plot((1:round(maxLag/2))*timeLapse,straightLine,'r--')

        %set axes limit
        axis([0 maxLag*timeLapse 0 1.2*max(timeAverage(:,1))]);
        
        %write axes labels
        xlabel('Lag (s)');
        ylabel('Time average (\mum^2)');
        title('From original coordinates')

        %write legend
        text(maxLag*timeLapse/2,0.4*max(timeAverage(:,1)),...
            sprintf(' Black: Kinetochore motion \n Red: Pure Brownian motion'));

        %hold off figure
        hold off

    end
    
end

%% ~~~ the end ~~~ %%

