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
potentialSlaves = {'forceMag','amp2','amp3'};
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
persistent imgStack tMap imgStack2 tracksNA curChanPath
%% Load tracksNA
finalProc = faPackage.getProcess(11);

%% Load the idsClassifiedStruct
iClasses = finalProc.loadChannelOutput(iChan,'output','idClass');
idGroupLabel= 1*iClasses.idGroup1 + ...
                2*iClasses.idGroup2 + ...
                3*iClasses.idGroup3 + ...
                4*iClasses.idGroup4 + ...
                5*iClasses.idGroup5 + ...
                6*iClasses.idGroup6 + ...
                7*iClasses.idGroup7 + ...
                8*iClasses.idGroup8 + 9*iClasses.idGroup9;

%% Load imgStack, forceStack and anyother stack if it exists.
if isempty(curChanPath) || ~strcmp(curChanPath, MD.channels_(1).channelPath_)
    curChanPath = MD.channels_(1).channelPath_;
    tracksNA=finalProc.loadChannelOutput(iChan,'output','tracksNA');
    [imgStack, tMap, imgStack2] = getAnyStacks(MD);
end

%% Launch pickAdhesion window with labeled adhesions with a right color and
% unlabed ones with white color. Get the right classes per newly selected
% adhesions
if PickManually
    waitHan = msgbox({'You will see the cell window and classified adhesion tracks with color';
    'label. After closing this window, you can select colored (classified)';
    'adhesions to see their time-series and associated assembly rates'});

    uiwait(waitHan);    
    pickAdhesionTracksInteractive(tracksNA, imgStack,...
        'movieData',MD,'tMap',tMap, 'imgMap2',imgStack2, 'idSelected',iClasses);
end
%% Representative case(s)
if RepClass==0
    RepClass=1:9;
end
numAvgWindow=0;
preDetecPeriod = 60; %seconds

for curClass=RepClass
    % Get the tracks of the same class
    curClassTracks = tracksNA(idGroupLabel==curClass);
    % Get the distribution of the time lag
    curFirstIncreseTimeIntAgainstSlave = ...
        calculateFirstIncreaseTimeTracks(curClassTracks,numAvgWindow,...
        preDetecPeriod,tInterval,'slaveSource',potentialSlaves{iSlave});
    
    % Get the median
    medianLag = nanmedian(curFirstIncreseTimeIntAgainstSlave);
    disp(['The median of the time lags of ' potentialSlaves{iSlave} ' of the class ' num2str(curClass) ...
        ' against amp1 is ' num2str(-medianLag,2) ' sec.'])
    disp(['The current movie is: ' MD.getFullPath '.'])
    % Get the repTracks where their time lags are the same as the median or
    % median + shift.
    trackID = curFirstIncreseTimeIntAgainstSlave == medianLag - TimeShiftFromMedian;
    % because the function generates a number of the master against the
    % slave (but usually we want the slave against the master), we subtract
    % (instead of add) the timeshift.
    
    % Loop through each id and generate the summary figures
    gPath = [faPackage.outputDirectory_ filesep 'RepTracks_Class' num2str(curClass)];
    if ~exist(gPath,'dir')
        mkdir(gPath)
    end
    for ii=find(trackID')
%         try
            h2 = showSingleAdhesionTrackSummary(MD,curClassTracks(ii),imgStack,tMap,imgStack2, ii,gPath);
            close(h2)
%         catch
%             disp(['error on this track: ' num2str(ii)])
%             continue
%         end
    end

end
disp(['The number of tracks identified: ' num2str(sum(trackID))])
disp(['Figures will be generated and stored in FocalAdhesionPackage/RepTracks_Class' num2str(curClass) '.'])
    


