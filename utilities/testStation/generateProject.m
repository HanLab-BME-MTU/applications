function [synthMovie, slist, idlist, dataProperties,storeDir,data2FileName] = generateProject(projectNumber, positions, inputDataProperties)
%GENERATEPROJECT is a utility to generate a synthetic project directory
%
% INPUT projectNumber: number/name of movie (moviename: trackTest_#)
%       positions    : positions (and amplitudes) of tags;
%                       nTimepoints-by-4-by-nTags array. 
%                       Identical tag positions will lead to a fusion spot.
%       inputDataProperties : (opt) dataProperties structure with the
%                               fields you want to have changed.
%
% jonas 3/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% outline
% 1) check computer, and decide where the data could be found
% 2) check if the data is already there
% 3) if yes: load data, return
% 4) if no: 
%       - calculate slist, idlist from position array
%       - generate movie (noise-free)
%       - save all to file

%===============
% CHECK MACHINE
%===============

% we're on the laptop if we're in windows and "home" is D:
isLaptop = ispc && strcmp(getenv('HOME'),'G:\matlab');

% on the laptop: path is D:\matlab\trackTest
% on the workstation: path is {Biodata-1}\trackTest

if isLaptop
    dataPath = 'G:\matlab\trackTest';
else
    [bioDataDir, oldDir] = cdBiodata(0);
    % go to right dir
    cd(['..' filesep 'trackTest']);
    dataPath = pwd;
    cd(oldDir);
end

%===============

%==================
% TRY TO LOAD DATA
%==================

% if project number is not an integer, we are doing a variant. Load
% original data, then check if we need to copy stuff.
projectString = num2str(projectNumber);
pointIdx = findstr(projectString,'.');
if ~isempty(pointIdx)
    % keep original project string
    variantString = projectString(pointIdx:end);
    projectString = projectString(1:pointIdx-1);
else
    variantString = [];
end

try
    
    % try to cd to the directory. If not found, we will create a project
    storeDir = [dataPath filesep 'testMovie_' projectString];
    
    % check whether we have a variant. If yes and directory exists, go
    % there directly. If yes and no directory exists, we copy the whole
    % diretory.
    if ~isempty(variantString)
        storeDirV = [storeDir, variantString];
        if ~exist(storeDirV,'dir')
            copyfile(storeDir,storeDirV);
        end
        storeDir = storeDirV;
    end
            
    oldDir = cd(storeDir);
    
    
    % find data2-file
    listOfFiles = searchFiles('-data2-','log',pwd,0,'new');
    
    % load data2-file
    data2FileName = listOfFiles{1};
    load(data2FileName);
    
    % load idlist, dataProperties, slist, rawMovie
    load(data2File.synthIdlist); 
    load(data2File.synthSlist);
    load(data2File.dataProperties);
    synthMovie = cdLoadMovie({data2File.movieName,'synth'});
    
    cd(oldDir)
    
    % we're done here. return
    return
    
catch
    [errormsg, msgID] = lasterr;
    
    if strcmp(msgID,'MATLAB:cd:NonExistentDirectory') ||  ~isdir(storeDir)...
            || strcmp(msgID,'MATLAB:load:couldNotReadFile')
        % we continue
    else
        % something else went wrong
        rethrow(lasterror);
    end
    
end

%===================

%============================
% GENERATE DATA
%============================

% set correct movie size, transform  positions
sptmp = size(positions);
sp = ones(1,3);
sp(1:length(sptmp)) = sptmp;
inputDataProperties.movieSize(4) = sp(1);

% assign and possibly update dataProperties
dataProperties = defaultDataProperties(inputDataProperties);

% make a movie based on the positions (noise free)
synthMovie = quickSynthMovie(positions, dataProperties);
% make movie header from data properties
r3dMovieHeader = getMovieHeader(dataProperties);

% transform positions before creating idlist
pixelSize = repmat([dataProperties.PIXELSIZE_XY,...
    dataProperties.PIXELSIZE_XY,...
    dataProperties.PIXELSIZE_Z],[sp(1),1,sp(3)]);
positions(:,1:3,:) = positions(:,1:3,:) .* pixelSize;

% now use positions to create a "pseudo"-slist and an idlist
idlist = positions2Idlist(positions);
slist  = idlist2slist(idlist);



%==================================
% STORE DATA
%==================================

% Store data in -data2- format

% make new directory
oldDir = cd(dataPath);
projectName = ['testMovie_' projectString];
dataProperties.name = projectName;
storeDir = [dataPath,filesep,projectName];
if ~isdir(storeDir)
    mkdir(storeDir);
end
cd(storeDir);

% save slist, idlist, dataProperties, synthMovie, movieHeader
idName = ['idlist-', nowString];
slName = ['slist-', nowString];
movieName = [projectName,'.r3c'];

save(idName,'idlist');
save(slName,'slist');
writemat(movieName,synthMovie);
movieTimeString = nowString;
save('dataProperties','dataProperties');
save('r3dMovieHeader','r3dMovieHeader');

% save -data2- file
[data2File, data2FileName] = createData2File(projectName,pwd);

data2File.dataProperties = 'dataProperties';
data2File.r3dMovieHeader = 'r3dMovieHeader';
data2File.movieName      = movieName;
data2File.slist          = slName;
data2File.idlist         = idName;
data2File.synthSlist     = slName;
data2File.synthIdlist    = idName;
data2File.lastResult     = 'idlist';

% store data2File and append logfile
save(data2FileName,'data2File');
logID = fopen([data2FileName, '.log'],'w');
fprintf(logID, '%s: synthetic movie %s created\n',movieTimeString, movieName);
fprintf(logID, '%s: data2File updated with synthetic raw data\n',nowString);
fclose(logID);

% if variant: copy everything
if ~isempty(variantString)
    storeDirV = [storeDir, variantString];
    copyfile(storeDir,storeDirV);
    storeDir = storeDirV;
end

% go back to original directory
cd(oldDir);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idlist = positions2Idlist(positions)
%IDLISTFROMPOSITIONS generates an idlist from a positions-array

nTimepoints = size(positions,1);
nTags = size(positions,3);
idlist(1:nTimepoints) = struct('linklist',zeros(nTags,12),'centroid',zeros(1,3),'stats',[],'info',[]);
labelcolor = [repmat('tag_',[nTags,1]),num2str([1:nTags]')]; %string tag_  1 - tag_999
labelcolor = cellstr(labelcolor); % conversion here, because matlab7 regexprep does not work on string arrays
labelcolor = regexprep(labelcolor,' ',''); %labelcolor = tag_1 - tag_999
idlist(1).stats.labelcolor = labelcolor;
idlist(1).stats.labellist = labelcolor;
idlist(1).weight = repmat(0.5,nTimepoints,1);
idlist(1).stats.status = {};
idlist(1).stats.maxColor = 2^(nTags+1);
idlist(1).stats.created = date;
for t = 1:nTimepoints
    % fill amp and positions
    idlist(t).linklist(:,[10,9,11,8]) = squeeze(positions(t,:,:))';
    % write time
    idlist(t).linklist(:,1) = t;
    % write spotnumber
    [dummy,dummy,idlist(t).linklist(:,2)] = unique(idlist(t).linklist(:,9:11),'rows');
    % write tag color
    idlist(t).linklist(:,4) = 2.^[0:nTags-1]';
    % write spot color - it's about time we got rid of this
    for i = 1:nTags
        idlist(t).linklist(i,3) =...
            sum(idlist(t).linklist(idlist(t).linklist(:,2)==i,4));
    end
    % write linkup, linkdow
    if t>1
        idlist(t).linklist(:,6) = idlist(t-1).linklist(:,3);
        idlist(t-1).linklist(:,7) = idlist(t).linklist(:,3);
    end
    
    % do chi^2=1 and capture all the uncertainty in the q-matrix
    idlist(t).linklist(:,12)=1;
    
    % we can't write in the local noise variance, because there is no noise
    % in the movie yet!
    
    % centroid
    idlist(t).centroid = mean(idlist(t).linklist(:,9:11),1);
    
    % Q-matrices: put 1 here, change according to noise of movie
    qdiag = repmat(1,[1,3*nTags]);
    idlist(t).info.detectQ_Pix = diag(qdiag);
    idlist(t).info.trackQ_Pix = [];
    
end

%=========================================================================

function synthMovie = quickSynthMovie(positions, dataProperties)
%QUICKSYNTHMOVIE calculates a synthetic Gaussian movie

% collect info
nTags = size(positions,3);

% preassign movie
frameSize = dataProperties.movieSize(1:3);
nTimepoints = size(positions,1);
synthMovie = zeros([frameSize,1,nTimepoints]);

% calculate size of psf-patch - it has to be at least large enough to
% contain psf values of 1% of amplitude
psfSize = roundOddOrEven(...
    2*sqrt(-2*dataProperties.FILTERPRM(1:3).^2*log(0.01)),'odd','inf');
hPsf=floor(psfSize/2);

% Loop through timepoints, then through spots. Calculate PSF-patch (we need
% to because of sub-pixel shift), and put into frame
for t = 1:nTimepoints
    for i = 1:nTags
        
         % read position, transform to pix
            pos = positions(t,1:3,i);

            % calculate position within pixel
            posfullpix = floor(pos);
            
            % [0,0,0] is in the center of the center pixel
            subpixelShift = pos - posfullpix;
            
            % make gaussPatch and multiply by amplitude
             psfData = GaussMask3D(...
                 dataProperties.FILTERPRM(1:3),psfSize,subpixelShift) *...
             positions(t,4,i);
             
            % put into movie - add onto existing intensities
            synthMovie(:,:,:,:,t) = ...
                pushstamp3d(synthMovie(:,:,:,:,t),...
                psfData, posfullpix, 0, 1);
            
    end % loop tags
end % loop time


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r3dMovieHeader = getMovieHeader(dataProperties)

r3dMovieHeader.pixelX=dataProperties.PIXELSIZE_XY;
r3dMovieHeader.pixelY=dataProperties.PIXELSIZE_XY;
r3dMovieHeader.pixelZ=dataProperties.PIXELSIZE_Z;
r3dMovieHeader.lensID=dataProperties.LENSID;
r3dMovieHeader.numCols=dataProperties.movieSize(2);
r3dMovieHeader.numRows=dataProperties.movieSize(1);
r3dMovieHeader.numZSlices=dataProperties.movieSize(3);;
r3dMovieHeader.numTimepoints=dataProperties.movieSize(4);;
r3dMovieHeader.numWvs=1;
r3dMovieHeader.wvl=dataProperties.WVL;
time = dataProperties.frameTime';
r3dMovieHeader.Time=time(:)';
r3dMovieHeader.expTime=dataProperties.expTime;
r3dMovieHeader.ndFilter=dataProperties.NDfilter;
r3dMovieHeader.cropInfo=dataProperties.crop;
r3dMovieHeader.correctInfo=[];
