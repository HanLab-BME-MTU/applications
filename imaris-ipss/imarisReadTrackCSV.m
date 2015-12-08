function [trk,frameTimes] = imarisReadTrackCSV(positionFile,timeFile)
%%IMARISREADTRACKCSV reads tracking data from CSV files exported from imaris
%
% tracks = imarisReadTrackCSV
% [tracks,frameTimes] = imarisReadTrackCSV
% ... = imarisReadTrackCSV(positionFile)
% ... = imarisReadTrackCSV(positionFile,timeFile)
%
%
%Hunter Elliott
%10/2014

%% ---- Parameters ---- %%


%Field name conversion cells 
imsName = {'PositionX','PositionY','PositionZ','Time'};%Imaris name
outName = {'x',            'y',           'z','Frame'};%Output field name
nFieldCopy = numel(outName);
p.Verbose = true;%In case we want to make optional later.

%% -----  input ----- %%

ip = inputParser;
ip.addParamValue('Verbose',true,@islogical);%If true, text status/progress displayed

if nargin < 1
    positionFile = '';
end

[positionPath,positionFile] = optionalFileInput(positionFile,'*.csv','Select an imaris position CSV file:');

if nargin < 2
    timeFile = '';
end

if isempty(timeFile);
    %Try guessing the time file based on position file:
    timeFile = [positionPath strrep(positionFile,'Position','Time')];    
    
    if ~exist(timeFile,'file')
        %If guess failed, ask the user...
        timeFile ='';
    else
        if p.Verbose;disp(['Auto-found time file, using "' timeFile '"']);end
    end        
end

[timePath,timeFile] = optionalFileInput(timeFile,'*.csv','Select an imaris time CSV file:');

    
%% ----- File Reading ----- %%
%Read the raw track data using matlab auto-generated import sub-functions

warning('off','MATLAB:table:ModifiedVarnames');%Disable warning to prevent user panic...

if p.Verbose;disp('Reading CSV files...');end

posTab = readtable([positionPath positionFile],'Format','%f%f%f%s%s%s%f%f%f%[^\n\r]','HeaderLines',3);
timeTab = readtable([timePath timeFile],'Format','%f%s%s%f%f%f%[^\n\r]','HeaderLines',3) ;


%% ----- Track Splitting ---- %%
%Re-arrange the data into a more sensible format and separate individual
%tracks

%First get the frame/time correspondence
frameTimes = unique([timeTab.Time timeTab.Value],'rows');
%In case there were no tracks on some frames, we re-create this here...
dT = min(diff(frameTimes(:,2)));
maxFrame = max(frameTimes(:,1));
frameTimes = (0:dT:(dT*(maxFrame-1)))';


%Get all track IDs
allTrkID = unique(posTab.TrackID);
allTrkID = allTrkID(~isnan(allTrkID));%Remove any spots without tracks
nTrk = numel(allTrkID);

if p.Verbose;disp(['Found ' num2str(nTrk) ' tracks, separating....']);end


%TEMP - pre-initialize structure!

for j = 1:nTrk
    
    currRows = posTab.TrackID == allTrkID(j);
    
    if nnz(currRows) > 0
    
        %Copy position and frame info
        for k = 1:nFieldCopy                
            trk(j).(outName{k}) = posTab.(imsName{k})(currRows);                
        end        
        %Get actual time data    
        trk(j).Time = frameTimes(trk(j).Frame);
        %Specify units and ID    
        trk(j).Time_Units = timeTab.Unit{1};
        trk(j).Position_Units = posTab.Unit{1};
        trk(j).TrackID = allTrkID(j);
    end
    
    
end

if p.Verbose;disp(['Done!']);end

