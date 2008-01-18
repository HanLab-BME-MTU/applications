function polyDepolyTimeAvg(runInfo,nFrms2Avg,timeStepSize,startFrm,endFrm)
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
    startFrm=1; default
    disp(['polyDepolyTimeAvg: startFrm = ' num2str(startFrm)])
end

if nargin<5 || isempty(endFrm)
    endFrm=[]; % will assign to nTotalFrames 
    disp(['polyDepolyTimeAvg: endFrm will number of frames in directory'])
end

% create subdirectory for time averaging the spatially averaged maps
% (first remove old directory if it exists)
turnTmAvgDir=[runInfo.turnDir filesep 'turnTmAvgDir'];
if isdir(turnTmAvgDir)
    rmdir(turnTmAvgDir,'s');
end
tmMatsDir=[turnTmAvgDir filesep 'mapMats'];
mkdir(tmMatsDir);


% get total number of images for correct digits while naming
[listOfImages] = searchFiles('.tif',[],runInfo.imDir);
nImTot=size(listOfImages,1);
s=length(num2str(nImTot));
strg=sprintf('%%.%dd',s);
indxStr=sprintf(strg,nImTot);

% count frames for proper naming
turnSpAvgMatDir=[runInfo.turnDir filesep 'turnSpAvgDir' filesep 'mapMats'];
[listOfFiles] = searchFiles('polyDepoly',[],turnSpAvgMatDir);
nFrames=size(listOfFiles,1);
s=length(num2str(nFrames));
strg=sprintf('%%.%dd',s);

% figure out which frames to average in a given iteration
if isempty(endFrm)
    endFrm=nFrames;
end
firstFrm=[startFrm:timeStepSize:endFrm]; lastFrm=firstFrm+nFrms2Avg-1;
firstFrm(lastFrm>endFrm)=[];             lastFrm(lastFrm>endFrm)=[];

nIterations=length(firstFrm); % how many iterations needed to go thru all frame ranges
polyDepolySum=zeros(nIterations,2); % store [sum(Poly) sum(Depoly)] for each frame range
polyDepolyFrmStack=zeros(runInfo.imL,runInfo.imW,nFrms2Avg); % store the nFrms2Avg maps here

for i=1:nIterations
    counter=1;
    for j=firstFrm(i):lastFrm(i)
        indxStr=sprintf(strg,j);
        iMap=load([turnSpAvgMatDir filesep 'polyDepoly' indxStr '.mat']);
        polyDepolyFrmStack(:,:,counter)=iMap.polyDepoly;
        counter=counter+1;
    end
    
    % get average value of poly/depoly at each pixel over the set of frames
    polyDepolyTmAvg=nanmean(polyDepolyFrmStack,3); 
    indxStr1=sprintf(strg,firstFrm(i));
    indxStr2=sprintf(strg,lastFrm(i));
    save([tmMatsDir filesep 'polyDepolyTmAvg' indxStr1 '_' indxStr2 '.mat'],'polyDepolyTmAvg');

    polyDepolySum(i,1)=sum(polyDepolyFrmStack(polyDepolyFrmStack(:)>0)); % sum poly over set of frames
    polyDepolySum(i,2)=sum(polyDepolyFrmStack(polyDepolyFrmStack(:)<0)); % sum depoly over set of frames
end

% take the average integrated (summed) poly/depoly over the set of frames
polyDepolySumAvg=polyDepolySum./nFrms2Avg; 
indxStr1=sprintf(strg,startFrm);
indxStr2=sprintf(strg,nFrms2Avg);
indxStr3=sprintf(strg,endFrm);
save([tmMatsDir filesep 'tmAvgParams'],'polyDepolySumAvg','nFrms2Avg','timeStepSize','startFrm','endFrm');


