function [data,orientation] = calculateTrajectoryFromIdlist(idlist,dataProperties,tag1,tag2,opt)
%CALCULATETRAJECTORYFORMIDLIST calculates distance trajectories from an idlist
%
%SYNOPSIS [data, orientation] = calculateTrajectoryFromIdlist(idlist,dataProperties,tag1,tag2,mode)
%
%INPUT    idlist         : any type of idlist
%         dataProperties : the corresponding data properties
%         tag1,2         : number (1:nTag) of first and second tag, OR string from labelcolor
%   if the first four arguments are omitted, the user is asked for input
%         opt            : optional options structure 
%                           .info :   struct that should be returned in the
%                                     output field info. 
%                           .calc2d : [{0}/1/2] depending on whether the
%                                     normal 3D-data or the maxProjection
%                                     or the in-focus-slice data should be
%                                     used.
%
%OUTPUT   data           : structure containing trajectory data
%                           .time       = [time in sec, sigmaTime]
%                           .timePoints = corresponding timepoint in raw data movie
%                           .distance   = [distance in um (tag1-tag2), distanceSigma]
%                           .info: substructure with fields helping to
%                                  describe the data
%                               .tags   cell containing the two tags used
%                                       calculate the trajectory {'tag1';'tag2'}
%
%         orientation    : 1-by-n matrix of angles [rad] from the xy-plane
%
%
%REMARKS  fusions will not be considered as valid timepoints
%
%c: 11/03 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=========================
%-------TEST INPUT--------
%=========================

%=========================
%check nargin (we do not really check for correct idlist, dataProperties)
if nargin < 4 | isempty(idlist) | isempty(dataProperties) | isempty(tag1) | isempty(tag2)
    [fileName,pathName] = uigetfile({'*-data-??-???-????-??-??-??.mat','project data files'},'select project data file');
        
        if fileName==0;
            h = errordlg('No data loaded','Warning!');
            uiwait(h);
            return % end evaluation here
        end
        cd(pathName);
        data = load(fileName); % loads everything into the structure data
        if ~isfield(data,'dataProperties')
            h = errordlg('No dataProperties in project data: corrupt data file','Warning!');
            uiwait(h);
            return % end evaluation here
        else
            dataProperties = data.dataProperties;
        end
        %---let the user choose which idlist to load
        
        % find which idlists there are
        dataFieldNames = fieldnames(data);
        idnameListIdx = strmatch('idlist',dataFieldNames);
        idnameList = dataFieldNames(idnameListIdx);
        
        % have the user choose, if there is more than one entry left
        switch length(idnameList)
            case 0 %no idlist loaded. continue w/o loading
                
                idname = '';
                h = errordlg('No idlist found in project data','Warning!');
                uiwait(h);
                
            case 1 % only one idlist loaded. Continue
                
                idname = char(idnameList);
                
            otherwise %l et the user choose
                idSelect = chooseFileGUI(idnameList);
                if isempty(idSelect)
                    idname = '';
                    h = errordlg('No data loaded!','Warning!');
                    uiwait(h);
                    return % end evaluation here
                else
                    idname = idnameList{idSelect};
                end
        end
        
        eval(['idlist = data.',idname,';']);
        
        clear data
        
        % get tags
        tagList = idlist(1).stats.labelcolor;
        
        h = helpdlg('Please choose first tag','');
        uiwait(h);
        tagNum1 = chooseFileGUI(tagList);
        tag1 = tagList{tagNum1};
        tagList(tagNum1)=[];
        
        h = helpdlg('Please choose second tag','');
        uiwait(h);
        tagNum2 = chooseFileGUI(tagList);
        tag2 = tagList{tagNum2};
        
        clear tagNum1 tagNum2 tagList
end
%=========================

%check options structure
%set defaults
info = [];
calc2d = 0;

%go through fields of opt, change defaults with input values if they exist
if nargin < 5 | isempty(opt)
    % keep defaults
else
    if isfield(opt,'info')
        info = opt.info;
    end
    if isfield(opt,'calc2d')
        calc2d = opt.calc2d;
    end
end

%check tags and convert to numbers
try
    labelcolor = idlist(1).stats.labelcolor;
    if isstr(tag1)
        tag1 = find(strcmp(tag1,labelcolor));
    end
    if isstr(tag2)
        tag2 = find(strcmp(tag2,labelcolor));
    end
    if isempty(tag1)|isempty(tag2)|any([tag1;tag2]>length(labelcolor))
        error('no tags found')
    else
        tags = idlist(1).stats.labelcolor([tag1,tag2],1);
    end
    %define indexList for indexing Q-matrix
    tagIdxList = [(tag1-1)*3+1:tag1*3,(tag2-1)*3+1:tag2*3];
catch
    error(['bad idlist or pointer to nonexistent tag: ',lasterr])
end



%calc pix2mu-matrix
p2mM = diag([dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_Z]);
pix2muMat = blkdiag(p2mM,p2mM);

%=========================
%---END TEST INPUT--------
%=========================




%=========================
%-------CALC 2D-CASE IF SELECTED
%=========================
switch calc2d
    case 1 %maximum projection
        [outputData3S,outputDataMP] = threeVsTwoD(idlist, dataProperties);
        idlist = outputDataMP.idlist;
    case 2 %in-focus slices
        [outputData3S,outputDataMP] = threeVsTwoD(idlist, dataProperties);
        idlist = outputData3S.idlist;
    otherwise
        %keep normal idlist
end

if isempty(idlist)
    error('can not calculate trajectory: idlist is empty')
end
        
%=========================
%---END CALC 2D-CASE IF SELECTED
%=========================





%=========================
%-------READ COORDS AND SIGMAZERO AND Q
%=========================

%cat linklist and extract coords, chi2, timePoints
linkLists = cat(3,idlist.linklist);


%coords
coords = linkLists([tag1,tag2],9:11,:);

%already calc distanceVectors to find zero distances
distanceVectors = diff(coords,1,1);
%make it a n-by-3 vector
distanceVectors = squeeze(distanceVectors)';

%find fusions
fusIdx = find(sum(distanceVectors,2)==0);

%eliminate entries in distVec, linkLists
distanceVectors(fusIdx,:) = [];
linkLists(:,:,fusIdx) = [];

%chi (be backward compatible)
if size(linkLists,2)<12
    chi2 = [];
else
    chi2 = squeeze(linkLists([tag1,tag2],12,:))';
end

%timepoints
timePoints = squeeze(linkLists(1,1,:));


%loop through idlist to get Q-matrices (write 'anaDat' for sigma-calculation)
ct = 1;
anaDat(1:length(timePoints)) = struct('stats',[]);
for t = timePoints' %only look at good timepoints
    %get QMatrix and transform to microns
    if isfield(idlist(t).info,'trackQ_Pix') & ~isempty(idlist(t).info.trackQ_Pix) & isfield(idlist(t).info,'trackerMessage')
        sourceT = str2double(idlist(t).info.trackerMessage.source);
        anaDat(ct).stats.qMatrix = pix2muMat*(idlist(t).info.trackQ_Pix(tagIdxList,tagIdxList)+...
            idlist(sourceT).info.detectQ_Pix(tagIdxList,tagIdxList))*pix2muMat;  
    else
        anaDat(ct).stats.qMatrix = pix2muMat*idlist(t).info.detectQ_Pix(tagIdxList,tagIdxList)*pix2muMat;
    end
    if isempty(chi2)
        anaDat(ct).stats.noise = idlist(t).info.noise([tag1,tag2]);
    else
        anaDat(ct).stats.noise = chi2(ct,:);
    end
    ct = ct+1;
end

%=========================
%---END READ COORDS AND SIGMAZERO AND Q
%=========================



%=========================
%-------CALC DISTANCE & DISTANCESIGMA
%=========================

%calc norms
[distanceN,normedVectors]=normList(distanceVectors);

%calc sigma
distanceSigma = adgui_calcPlotData_distanceSigma(anaDat,[1 2],distanceVectors,distanceN,0);

%=========================
%---END CALC DISTANCE AND DISTANCESIGMA
%=========================



%=========================
%-------READ TIME AND TIMESIGMA
%=========================

timeAll = dataProperties.frameTime;
time = mean(timeAll(timePoints,:),2);
%timeSigma is half the time between (lastCol - firstCol of frameTime)
timeSigma = (timeAll(timePoints,end)-timeAll(timePoints,1))/2;

%=========================
%---END READ TIME AND TIMESIGMA
%=========================

% %-------READ CHI2 - don't need
% if isempty(chi2)
%     for i = 1:ct-1
%         chi2(i,:)=anaDat(i).stats.noise;
%     end
% end
% %---END READ CHI2



%=========================
if nargout > 1
    %orientation: [-pi/2...pi/2]. calculate from vN*[0 0 1]
    orientation = pi/2-acos(normedVectors(:,3));
else
    orientation = [];
end

%=========================
%-------ASSIGN OUTPUT
%=========================
data.distance   = [distanceN,distanceSigma];
data.time       = [time,timeSigma];
data.timePoints = timePoints;
data.info       = info;
data.info.tags  = tags;
% data.chi2 = mean(chi2,2);

%---END ASSIGN OUTPUT