function [] = showMostRepresentativeAssemblingAdhesion(MD,varargin)
%function [] = showMostRepresentativeAssemblingAdhesion(MD) searches
%through trajectories, pick one that shows bona fide new assembly and shows
%the best goodness of fit.
% input
%       MD                      MovieData file that has run FAPackage
%       RepClass                Specific class name (1-9). Type 0 for all classes and
%                               pick manually
%       PickManually            true if you want to pick manually (default: false)
%       TimeShiftFromMedian     time shift from the median time lag. e.g.
%                               The TimeLagFromMedian -1 is 1 sec the secondary signal (e.g. 
%                               vinculin) comes before the main signal (e.g. talin).
%       iSlave                  the slave channel. 1 is for TFM, 2 is for amp2, 3 is for amp3.
%                               The main is the one at channel 2 (e.g. talin-GFP channel) usually
% output
%       figures will be generated and stored in
%       FocalAdhesionPackage/RepTracks_Class
% Example: 
% showMostRepresentativeAssemblingAdhesion(MD,'RepClass',1,'TimeShiftFromMedian',-6, 'iSlave',2) 
% To sample adhesion track manually,
% showMostRepresentativeAssemblingAdhesion(MD,'RepClass',0,'PickManually',true,'iSlave',2)
% Sangyoon Han, March 4 2020

%% Input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('MD', @(x) isa(x,'MovieData'));
ip.addParameter('RepClass',2, @isnumeric); %0 for all classes and pick manually
ip.addParameter('PickManually',false, @islogical);
ip.addParameter('TimeShiftFromMedian',0, @isscalar); 
% e.g. The TimeLagFromMedian -1 is 1 sec the secondary signal (e.g.
% vinculin) comes before the main signal (e.g. talin).
ip.addParameter('iSlave',1, @isnumeric); 
% iSlave is the slave channel. 1 is for TFM, 2 is for amp2, 3 is for amp3.
% The main is the one at channel 2 (e.g. talin-GFP channel) usually.

ip.parse(MD,varargin{:});
RepClass=ip.Results.RepClass;
PickManually=ip.Results.PickManually;
TimeShiftFromMedian=ip.Results.TimeShiftFromMedian;
iSlave=ip.Results.iSlave;

tInterval = MD.timeInterval_;
potentialSlaves = {'forceMag','amp2','amp3','fret'};
% existingSlaveIDs = isfield(tracksNA,potentialSlaves);
if mod(TimeShiftFromMedian,tInterval)>0
    disp(['Previous time shift was ' num2str(TimeShiftFromMedian) '.'])
    TimeShiftFromMedian = tInterval*round(TimeShiftFromMedian/tInterval);
    disp(['The time shift was adjusted to ' num2str(TimeShiftFromMedian) ' to fit with the time interval.'])
end
%% Load FA package
faPackage=MD.getPackage(MD.getPackageIndex('FocalAdhesionPackage'));
% Load classification process
classProc = faPackage.getProcess(8);
iChan = find(classProc.checkChannelOutput);

%% persistent set up for large memory-requiring variable
persistent imgStack tMap imgStack2 fretMap tracksNA curChanPath
%% Load tracksNA
finalProc = faPackage.getProcess(11);

%% Load the idsClassifiedStruct
% iClasses = finalProc.loadChannelOutput(iChan,'output','idClass');
iClassObj = load([finalProc.funParams_.OutputDirectory filesep 'data' ...
    filesep 'idGroups.mat'],'idGroups');
iClasses = iClassObj.idGroups;
idGroupLabel = 1*iClasses{1};
for ii=2:9 
    idGroupLabel = idGroupLabel+ ii*iClasses{ii};
end
% idGroupLabel= 1*iClasses.idGroup1 + ...
%                 2*iClasses.idGroup2 + ...
%                 3*iClasses.idGroup3 + ...
%                 4*iClasses.idGroup4 + ...
%                 5*iClasses.idGroup5 + ...
%                 6*iClasses.idGroup6 + ...
%                 7*iClasses.idGroup7 + ...
%                 8*iClasses.idGroup8 + 9*iClasses.idGroup9;

%% Load imgStack, forceStack and anyother stack if it exists.
if isempty(curChanPath) || ~strcmp(curChanPath, MD.channels_(1).channelPath_) ...
    || numel(tracksNA) ~= sum(idGroupLabel==RepClass)
    if RepClass>0
        tracksNA=finalProc.loadChannelOutput(iChan,'output','tracksNA','idSelected',find(idGroupLabel==RepClass)');
    else
        tracksNA=finalProc.loadChannelOutput(iChan,'output','tracksNA');
    end
end
%% Representative case(s)
if RepClass==0
    iRepClass=1:9;
else
    iRepClass = RepClass;
end
numAvgWindow=0;
preDetecPeriod = 60; %seconds

% for curClass=iRepClass
curClass=iRepClass;
    % Get the tracks of the same class
%     curClassTracks = tracksNA; %(idGroupLabel==curClass);
    % Read the intensity again
%     tracksNA=readIntensityFromTracks(tracksNA,imgStack,1,'extraLength',120,'movieData',MD,'reTrack',true);

    % Get the distribution of the time lag
    if iSlave>0
        curFirstIncreseTimeIntAgainstSlave = ...
            calculateFirstIncreaseTimeTracks(tracksNA,numAvgWindow,...
            preDetecPeriod,tInterval,'slaveSource',potentialSlaves{iSlave});

        % Get the median
        medianLag = nanmedian(curFirstIncreseTimeIntAgainstSlave);
        disp(['The median of the time lags of ' potentialSlaves{iSlave} ' of the class ' num2str(curClass) ...
            ' against amp1 is ' num2str(-medianLag,2) ' sec.'])
    end
    disp(['The current movie is: ' MD.getFullPath '.'])
    % Get the repTracks where their time lags are the same as the median or
    % median + shift.
    if iSlave>0
        trackID = curFirstIncreseTimeIntAgainstSlave == medianLag - TimeShiftFromMedian;
        % because the function generates a number of the master against the
        % slave (but usually we want the slave against the master), we subtract
        % (instead of add) the timeshift.

        % Exporting
        lifetimes = arrayfun(@(x) x.lifeTime, tracksNA);
        t = table(curFirstIncreseTimeIntAgainstSlave,lifetimes);
        writetable(t,[finalProc.funParams_.OutputDirectory filesep 'data' filesep 'timeLagG' num2str(RepClass) '.xlsx']);
        save([finalProc.funParams_.OutputDirectory filesep 'data' filesep 'timeLagG' num2str(RepClass) '.mat'],...
            'curFirstIncreseTimeIntAgainstSlave','lifetimes');
    end    
    
    if isempty(imgStack) || ~strcmp(curChanPath, MD.channels_(1).channelPath_)
        [imgStack, tMap, imgStack2, fretMap] = getAnyStacks(MD);
        curChanPath = MD.channels_(1).channelPath_;
    end
    %% Launch pickAdhesion window with labeled adhesions with a right color and
    % unlabed ones with white color. Get the right classes per newly selected
    % adhesions
    if PickManually
        waitHan = msgbox({'You will see the cell window and classified adhesion tracks with color';
        'label. After closing this window, you can select ';
        ', by using Data Tips icon, a colored (classified)';
        'adhesion to see their time-series and associated assembly rates'});

        uiwait(waitHan);    
        pickAdhesionTracksInteractive(tracksNA, imgStack,...
            'movieData',MD,'tMap',tMap, 'imgMap2',imgStack2,'fretMap',fretMap, 'idSelected',iClasses);
    end

    %% background substraction
    if ~isempty(imgStack)
        imgClass = class(imgStack);
        imgStackBS=zeros(size(imgStack));
        for ii=1:size(imgStack,3)
            curImg=imgStack(:,:,ii);
            imageBackground = filterGauss2D(curImg,30);
            %calculate noise-filtered and background-subtracted image
            imgStackBS(:,:,ii) = curImg - cast(imageBackground,imgClass);
        end
    end
    if ~isempty(imgStack2)
        imgClass = class(imgStack2);
        imgStackBS2=zeros(size(imgStack2));
        for ii=1:size(imgStack2,3)
            curImg=imgStack2(:,:,ii);
            imageBackground = filterGauss2D(curImg,30);
            %calculate noise-filtered and background-subtracted image
            imgStackBS2(:,:,ii) = curImg - cast(imageBackground,imgClass);
        end
    else
        imgStackBS2 = [];
    end

    
    % Loop through each id and generate the summary figures
    gPath = [faPackage.outputDirectory_ filesep 'RepTracks_Class' num2str(curClass)];
    if ~exist(gPath,'dir')
        mkdir(gPath)
    end
    if RepClass>1
        for ii=1:numel(tracksNA)
            h2 = showSingleAdhesionTrackSummary(MD,tracksNA(ii),imgStack,tMap,imgStack2, ii,gPath,[],imgStackBS,imgStackBS2);
            close(h2)
        end
    else
        if exist('trackID','var')
            for ii=find(trackID')
                h2 = showSingleAdhesionTrackSummary(MD,tracksNA(ii),imgStack,tMap,imgStack2, ii,gPath,[],imgStackBS,imgStackBS2);
                close(h2)
            end
        end
    end
% end

if exist('trackID','var')
    disp(['The number of tracks identified: ' num2str(sum(trackID))])
end
disp(['Figures will be generated and stored in FocalAdhesionPackage/RepTracks_Class' num2str(curClass) '.'])
    


