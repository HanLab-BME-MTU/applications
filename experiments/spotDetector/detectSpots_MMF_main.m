function [slist, debugData] = detectSpots_MMF_main(rawMovie,cord,dataProperties,testRatios,verbose,debug)
%DETECTSPOTS_MMF_MAIN is the main function of the mixture model fitting
%
% SYNOPSIS: [slist, optData] = detectSpots_MMF_main(rawMovie,cord,dataProperties,testRatio,verbose,options)
%
% INPUT rawMovie: raw movie frames
%       cord: coordinates array
%		dataProperties: data properties structure
%		testRatio: output of detectSpots_MMF_findAmplitudeCutoff. Can be
%                  empty
%		verbose: 0: no window, 1: waitbar, 2 or string: progressText
%		debug: debug status
%
% OUTPUT slist: nTimepoints-by-1 structure with fields:
%			.sp spots structure with fields...
%			debugData: optional additional output
%
% REMARKS
%
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: Jonas Dorn
% DATE: 07-Feb-2006
%
%

%=========================
% CHECK INPUT
%=========================

if nargin < 5 || isempty(verbose)
    verbose = 1;
end
if ischar(verbose)
    progTxt = verbose;
    verbose = 2;
else
    progTxt = '';
end
if nargin < 6 || isempty(debug)
    debug = 0;
end

% preassign debugData
debugData = [];

%=========================


%=========================
% PREPARE FITTING
%=========================

% find length of movie
nTimepoints = size(rawMovie,5);
frameSize = dataProperties.movieSize(1:3);

% find goodTimes, which timepoints analyzed if possible
if isempty(testRatios)
    % assign goodTimes, trueTime with what we know
    goodTimes = (1:nTimepoints)';
    trueTime = [1,nTimepoints];

else
    % find good times from testRatio
    goodTimes = find(~cellfun('isempty',testRatios));

    if ~isempty(goodTimes)
        % find which timepoints we're actually analyzing
        tmp = cat(1,testRatios{:});
        trueTime(2) = max(tmp(:,1));
        trueTime(1) = min(tmp(:,1));
    else
        % set trueTime to -99 and  continue. There will be an empty slist
        trueTime = [-99,-99];
    end
end



% preassign slist
slist(1:nTimepoints) = struct('sp',[],'statistics',[],'parms',[],...
    'COM',{cord.COM});

% turn off singularMatrix warning
warningState = warning;
warning off MATLAB:nearlySingularMatrix


switch verbose
    case 1
        h= mywaitbar(0,[],nTimepoints,...
            sprintf('Mixture Model Fitting, frames %i:%i',...
            trueTime(1),trueTime(2)));
    case 2
        progressText(0,progTxt);
end

%==========================


%==============================
% MIXTURE MODEL FITTING LOOP
%==============================

% init debug
% 1: collect fStats
if any(debug) == 1
    debugData.fStats = cell(nTimepoints,1);
end


% loop through time
for t=goodTimes'
    Q=[];
    qAmp = [];
    chi=[];
    snr=[];
    ct=1;

    % we need to keep the if clause here in case there is no testRatios
    if isfield(cord(t),'sp') &&  ~isempty(cord(t).sp)
        cordList=cat(1,cord(t).sp.cord);
        %change to matlab coords
        tc2=cordList(:,2);
        cordList(:,2)=cordList(:,1);
        cordList(:,1)=tc2;
        currentFrame=rawMovie(:,:,:,1,t);

        % DEBUG
        % nspots=size(cordList,1);
        % disp(sprintf('t : %i nsp : %i', t,nspots));

        % remove spots that don't pass the first round of fitting, anyway
        if ~isempty(testRatios)
            currentTest = testRatios{t};
            goodIdx = currentTest(:,2) >= dataProperties.amplitudeCutoff;
        else
            goodIdx = 1:size(cordList,1);
        end
        % if everything has been removed, it is bad. Therefore, we put at
        % least one of the coordinates back into the list. Fitting one spot
        % gives us an estimate for the background and ensures that all
        % fields will be filled. Also, just fitting one single spot does
        % not take up that much time.
        if ~any(goodIdx)
            goodIdx = 1;
        end
        cordList = cordList(goodIdx,:);

        while (~isempty(cordList))

            % find (again) spots that have to be fitted jointly
            [spotsidx, mask] = discernspots(cordList,frameSize,dataProperties);

            idxList=find(mask);
            mskData=currentFrame(idxList);

            % do the mixture-model fitting
            [numDist,ncordList,ampList,bg,statistics, dbTmp]=...
                fitTest(mskData,cordList(spotsidx,:),idxList,...
                frameSize,dataProperties, debug);
            
            % !!! HACK
            % somehow for Ipl1-bilobe-data, fitTest places the spots way
            % off. I have no idea why and no time to find out. Thus, use
            % cordList if negative amplitudeCutoff
            if dataProperties.amplitudeCutoff < 0
                ncordList = cordList(spotsidx,:);
            end
            
            % DEBUG - STORE TESTSTAT IN BASE
%             testStat = cat(2,dbTmp.testValue)';
%             testStat = [ones(length(testStat),1)*t,testStat];
%             ts = getappdata(0,'testStat');
%             ts = [ts;testStat];
%             setappdata(0,'testStat',ts)

            % read debugData
            if debug == 1
                debugData.fStats{t} = ...
                    [debugData.fStats{t};cat(1,dbTmp.fStats)];
            end

            if isempty(ncordList)
                statistics.parms=[];
                statistics.multi=0;
                statistics.chi=[];
                statistics.Qxx=[];
                statistics.snr=[];
                statistics.qAmp=[];
            end;
            %define new cordLst
            tempidx=1:size(cordList,1);
            tempidx(spotsidx)=0;
            tempidx=nonzeros(tempidx);
            cordList=cordList(tempidx,:);

            % add Q of spot(s)
            Q=blkdiag(Q,statistics.Qxx);
            qAmp = blkdiag(qAmp,statistics.qAmp);
            snr=[snr statistics.snr];

            %insert coordinates

            for i= 1:numDist
                % insert new center coords
                patch=stamp3d(currentFrame,[2 2 2],round(ncordList(i,:)));
                %change back to standard coords
                slist(t).sp(ct).cord=[ncordList(i,2) ncordList(i,1) ncordList(i,3)];
                slist(t).sp(ct).amp=ampList(i);
                slist(t).sp(ct).bg=bg;
                slist(t).sp(ct).mnint=mean(patch(:));
                % overlapping or single spot
                slist(t).sp(ct).mult=statistics.multi;
                ct=ct+1;
                chi=[chi statistics.chi];
            end; % for i= 1:numDist
        end; %  while (~isempty(cordList))


        % test for zero distance
        slist(t).parms=statistics.parms;

        if ~isempty(Q)
            % normally, there should not be a problem here anymore - the
            % tags that were close were fitted jointly and tested already.
            % However, there is the unlikely event that we had single tags
            % that were really badly fitted and thus not significantly
            % different from each other
            slist(t) = testdistance(slist(t),Q,chi,snr,dataProperties);
            % just add qAmp for the moment; don't check for removed
            % spots
            slist(t).statistics.qAmp = qAmp;
        else
            % the if and else here  is only because the nspots>=MAXSPOTS
            slist(t).parms=[];
        end;

    end;
    % progressbar
    switch verbose
        case 1
            mywaitbar(t/nTimepoints,h,nTimepoints);
        case 2
            progressText(t/nTimepoints);
    end
end; % loop time

%========================

%========================
% CLEAN UP
%========================

warning(warningState)

if verbose==1
    close(h);
end

%========================