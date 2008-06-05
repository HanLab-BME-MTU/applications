function polyDepolyTimeAvg(runInfo,nFrms2Avg,timeStepSize,startFrm,endFrm,movieName,mkMov)
%POLYDEPOLYTIMEAVG: time average kinetic maps and calculate mean total poly/depoly over frame range
%
% SYNOPSIS: [polyDepolySumAvg]=polyDepolyTimeAvg(runInfo,nFrms2Avg,timeStepSize,startFrm,endFrm)
%
% INPUT: runInfo        : structure containing a field with the turnover
%                         analysis directory path (from polyDepolyMap)
%        nFrms2Avg      : frame interval to average
%        timeStepSize   : number of frames between averages
%        startFrm       : first frame to use in calculation
%        endFrm         : last frame to use in calculation
%
%        e.g.) if you want frames 1-5 averaged, 3-8, 5-11, etc. then
%                         nFrms2Avg=5 and timeStepSize=3
%
% OUTPUT: /analysis/turn/turnTmAvgDir/mapMats : directory with .mat files
%              containing the averaged poly/depoly maps
%         polyDepolySumAvg : n x 2 matrix containing the average integrated
%              value for [poly depoly] calculated from raw,
%              spatially-averaged (over 4 grids) maps. The idea
%              is to compare the average amount of poly vs.
%              depoly over parts of the movie. This calculation
%              is not normalized for the number of pixels used
%              to make the measurement in each frame.  This value gets
%              saved in the mapMats folder in tmAvgParams.
%
% USERNAME: kathomps
% DATE: 11-Jan-2008

if nargin<1 || ~isstruct(runInfo) || ~isfield(runInfo,'turnDir')
    error('polyDepolyTimeAvg: runInfo should be a structure with field turnDir')
end

if nargin<2 || isempty(nFrms2Avg)
    nFrms2Avg=5; % default
    disp(['polyDepolyTimeAvg: nFrms2Avg = ' num2str(nFrms2Avg)])
end

if nargin<3 || isempty(timeStepSize)
    timeStep=1; % default
    disp(['polyDepolyTimeAvg: timeStep = ' num2str(timeStepSize)])
end

if nargin<4 || isempty(startFrm)
    startFrm=1; % default
end

if nargin<5 || isempty(endFrm)
    endFrm=[]; % will assign to nTotalFrames
end

% create subdirectory for time averaging the spatially averaged maps
% (first remove old directory if it exists)
turnTmAvgDir=[runInfo.turnDir filesep 'turnTmAvgDir_' num2str(nFrms2Avg) '_' num2str(timeStepSize)];
if isdir(turnTmAvgDir)
    rmdir(turnTmAvgDir,'s');
end
tmMatsDir=[turnTmAvgDir filesep 'mapMats'];
mkdir(tmMatsDir);
% make subdirectory for storing max-edge cell masks over the time interval
tmCellMasksDir=[tmMatsDir filesep 'cell_mask'];
mkdir(tmCellMasksDir);

% create subdirectory for histograms showing the kinetic activity from the
% time-averaged data
% (first remove old directory if it exists)
histDir=[turnTmAvgDir filesep 'activityHist'];
mkdir(histDir);

% count frames for proper naming
turnMapDir=[runInfo.turnDir filesep 'polyDepolyMaps'];
[listOfFiles] = searchFiles('polyDepoly',[],turnMapDir);
nFrames=size(listOfFiles,1);
s=length(num2str(nFrames));
strg=sprintf('%%.%dd',s);

% if 0, than average over the whole movie
if nFrms2Avg==0
    nFrms2Avg=nFrames;
    endFrm=nFrames;
end

% figure out which frames to average in a given iteration
if isempty(endFrm)
    endFrm=nFrames;
end
firstFrm=[startFrm:timeStepSize:endFrm]; lastFrm=firstFrm+nFrms2Avg-1;
firstFrm(lastFrm>endFrm)=[];             lastFrm(lastFrm>endFrm)=[];

nIterations=length(firstFrm); % how many iterations needed to go thru all frame ranges
avgPolyDepolyPerFrame=zeros(nIterations,2); % store [sum(avgPoly) sum(avgDepoly)] for each frame range
polyDepolyFrmStack=zeros(runInfo.imL,runInfo.imW,nFrms2Avg); % store the nFrms2Avg maps here

% use the cell masks of the images to average
cmDir=[runInfo.anDir filesep 'edge' filesep 'cell_mask'];
[listOfCellMasks] = searchFiles('.tif',[],cmDir,0);
roiMask=runInfo.roiMask; % max cell proj bounded by fieldGeom

if mkMov==1
    wholeMovieName=[histDir filesep 'actHist_' movieName '_' num2str(nFrms2Avg) '.avi'];
    aviobj = avifile(wholeMovieName,'fps',10);
end

for i=1:nIterations
    counter=1;
    cellMaskTmAvg=zeros(runInfo.imL,runInfo.imW);
    for j=firstFrm(i):lastFrm(i)
        % put the polyDepoly maps for the interval into a big array
        indxStr=sprintf(strg,j);
        iMap=load([turnMapDir filesep 'polyDepoly' indxStr '.mat']);
        polyDepolyFrmStack(:,:,counter)=iMap.polyDepoly;
        counter=counter+1;

        % get max edge projection cell mask
        cMask=double(imread([char(listOfCellMasks(j,2)) filesep char(listOfCellMasks(j,1))]));
        cellMaskTmAvg = cellMaskTmAvg | cMask;
    end
    % bound cell mask by fieldGeom
    cellMaskTmAvg=logical(cellMaskTmAvg.*roiMask);

    % get average value of poly/depoly at each pixel over the set of frames
    polyDepolyTmAvg=nanmean(polyDepolyFrmStack,3);
    indxStr1=sprintf(strg,firstFrm(i));
    indxStr2=sprintf(strg,lastFrm(i));
    save([tmMatsDir filesep 'polyDepolyTmAvg' indxStr1 '_' indxStr2 '.mat'],'polyDepolyTmAvg');

    % write the tif image of the max projection cell mask
    imwrite(cellMaskTmAvg,[tmCellMasksDir filesep 'timeAvgCellMask' indxStr1 '_' indxStr2 '.tif']);

    cellMaskTmAvg=swapMaskValues(cellMaskTmAvg,0,NaN); % 1's inside max cell; nan's outside
    maskedAvgActivity=polyDepolyTmAvg.*cellMaskTmAvg; % mask out the time averaged activity map
    avgPolyDepolyPerFrame(i,1)=sum(maskedAvgActivity(maskedAvgActivity(:)>0)); % sum of the time-averaged poly map
    avgPolyDepolyPerFrame(i,2)=sum(maskedAvgActivity(maskedAvgActivity(:)<0)); % sum of the time-averaged depoly map

    % create histogram of poly/depoly values from the set of nFrms2Avg
    temp=maskedAvgActivity(:);
    temp(isnan(temp))=[];
    hist(temp,100);
    axis([-0.04 0.04 0 2000]);
    saveas(gca,[histDir filesep 'actHist' indxStr1 '_' indxStr2 '.png']);
    frame = getframe(gca);
    if mkMov==1
    aviobj = addframe(aviobj,frame);
    end
end
if mkMov==1
    aviobj = close(aviobj);
end
save([tmMatsDir filesep 'tmAvgParams'],'avgPolyDepolyPerFrame','nFrms2Avg','timeStepSize','startFrm','endFrm');



