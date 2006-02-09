function [cordOut, debugData]=findoverlap(data,cord,dataProperties,verbose)
%FINDOVERLAP locates fluorescent tags in 3D data
%
% SYNOPSIS cords=findoverlap(data,cord)
%
% INPUT img   : stack time series
%       verbose: (opt) If 1 (default), waitbar is displayed.
%
% OUTPUT cord : center coordinates

% c: 18/6/01	dT

% optional input
if nargin < 4 || isempty(verbose)
    verbose = 1;
end

debugData = [];

%CONST DEFINITIONS
MAXSPOTS=dataProperties.MAXSPOTS;

% turn off singularMatrix warning
warningState = warning;
warning off MATLAB:nearlySingularMatrix

maxNumSpots = 99999; % maximum number of spots that can be fitted simultaneously
maxNumCoords = 99999; % maximum number of coordinates that can be checked by discernspots
tsteps=size(data,5);
if verbose
    h= mywaitbar(0,[],tsteps*2,'Checking for multiple spots...');
end

% init cordOut. Code changed to no longer overwrite input cord, but to
% build a new structure from scratch. This could potentially be a problem
% with huge structures, but it is much less error-prone
cordOut(1:tsteps) = struct('sp',[],...
    'statistics',[],'parms',[],'COM',{cord.COM});

%===================
% PREPARE FITTING
%===================

% find cutoff for amplitude. Make it the easy way first (we can try to
% program it faster later): For all frames, fit the local maxima (just N
% spots) and store amp/sigmaAmp. Make a histogram of that to find the best
% cut-off for good spots vs. noise. Then do the real fitting.
% Actually, we can already remove some bad fits, although we'll still
% better go through the entire fitTest for the moment.

% define store variable for testRatios
testRatios = cell(tsteps,1);

% loop t to fit local maxima
for t = 1:tsteps
    % find spots to fit
    if isfield(cord(t),'sp') &  ~isempty(cord(t).sp)
        %
        % have try-catch loop to watch out for memory issues
        try

            cordList=cat(1,cord(t).sp.cord);
            %change to matlab coords
            tc2=cordList(:,2);
            cordList(:,2)=cordList(:,1);
            cordList(:,1)=tc2;
            
            imgStk=data(:,:,:,1,t);
            
            % preassign current testRatios
            nSpots = size(cordList,1);
            testRatios{t} = zeros(nSpots,2);
            spotIdxList = 1:nSpots;
            
            while (~isempty(cordList))
                discerned = 0;
                if length(cordList) >= maxNumCoords
                    error('MATLAB:nomem',...
                        'too many coords (%i) - Matlab will run out of memory!',...
                        length(cordList));
                end
                [spotsidx mask] = ...
                    discernspots(cordList,size(imgStk),dataProperties);
                currentSpots = spotIdxList(spotsidx);
                
                discerned = 1;
                if length(spotsidx) >= maxNumSpots
                    error('MATLAB:nomem',...
                        'too many spots (%i) - Matlab will run out of memory!',...
                        length(spotsidx));
                end
                idxList=find(mask);
                mskData=imgStk(idxList);

                % do the maximum fitting
                testValue =...
                    fitTestFirstRound(mskData,cordList(spotsidx,:),idxList,...
                    size(imgStk),dataProperties);
                % collect testValues. Be careful:
                % - discernspots picks groups that can be any combination
                % - fitTest reverses the order of the spots
                currentRatios = testRatios{t};
                currentRatios(currentSpots(end:-1:1),:) = ...
                    [ones(length(testValue),1)*t,testValue(:)];
                testRatios{t} = currentRatios;
                
                % remove spots from cordList
                cordList(spotsidx,:)=[];
                spotIdxList(spotsidx) = [];

            end
        catch
            % check the type of error
            err = lasterror;
            if strmatch(err.identifier,'MATLAB:nomem')
                % we have a memory problem
                % remember length of spotsIdx
                if discerned
                    maxNumSpots = min(maxNumSpots,length(spotsidx));
                else
                    maxNumCoords = min(maxNumCoords,length(cordList));
                end
                
                % remove the spots that have already been analyzed
                testRatios{t} = [];

                disp(sprintf('MMF aborted in frame %i:\n%s',t,err.message))
            else
                % this is not a memory problem. rethrow the error
                rethrow(err)
            end
        end
    end % if ~isempty
    if verbose
        mywaitbar(t/tsteps/2,h,tsteps*2);
    end
end % loop time

% get cutoff. Transform cell to matrix and read out ratios only
tmp = cat(1,testRatios{:});
goodTimes = unique(tmp(:,1));
ratios = tmp(:,2);
clear tmp

% remove very small values from ratios. If there are multiple almost zeros,
% for example, there will be a peak at 0 that will screw up
% cutFirstHistMode.
ratios(ratios < 1e-3) = [];

% find cutoff

% first guess via cutFirstHistMode
[cutIdx, cutVal,sp] = cutFirstHistMode(ratios,0);

% now check the local minima in the vicinity of the cutoff
spder = fnder(sp);
zeroList = fnzeros(spder);
zeroList = zeroList(1,:);
% evaluate
zeroVals = fnval(sp,zeroList);

% look in zeroList. Find one value before cutVal, three after. Go into
% zeroVals and find lowest minimum
[dummy,closestIdx] = min(abs(zeroList - cutVal));

% check only the minimas that are close by; two to the right and
% one to the left (don't forget that between two minima there will
% always be a maximum!)
indexList = (closestIdx-2):(closestIdx + 4);
indexList(indexList < 1 | indexList > length(zeroVals)) = [];

% find lowest
[dummy, cutIdx] = min(zeroVals(indexList));
% and determine break value
cutValue = zeroList(indexList(cutIdx));

% plot cutoff
if 1
    figure('Name','cutoff'),
    axesH(1)=subplot(1,2,1);
    plot(ratios,'+'),
    axesH(2)=subplot(1,2,2);
    cutFirstHistMode(axesH(2),ratios);
    set(axesH(1),'NextPlot','add');
    plot(axesH(1),[1,length(ratios)],[cutVal,cutVal],'g')
    plot(axesH(1),[1,length(ratios)],[cutValue,cutValue],'r')
    set(axesH(2),'NextPlot','add');
    plot(axesH(2),[cutVal,cutVal],[0,100],'g',...
        [cutValue,cutValue],[0,100],'r')
end

% store cutValue in dataProperties
dataProperties.amplitudeCutoff = cutValue;


%===================================
% MAIN FITTING LOOP
%===================================

for t=goodTimes'
    Q=[];
    qAmp = [];
    chi=[];
    snr=[];
    ct=1;
    

    %======DEBUG

    %             if any(t == [8])
    %                 disp('press dbstep to continue!')
    %                 keyboard
    %             end
    %======DEBUG

    %     % if there are too many local maxima, the code could crash
    %     % due to memory issues. therefore, try/catch to have only one bad
    %     % frame, not lose an entire movie
    %
    %     try
    % if no spots in this time frame -> go back
    if isfield(cord(t),'sp') &  ~isempty(cord(t).sp)
        cordList=cat(1,cord(t).sp.cord);
        %change to matlab coords
        tc2=cordList(:,2);
        cordList(:,2)=cordList(:,1);
        cordList(:,1)=tc2;
        nspots=size(cordList,1);
        %disp(sprintf('t : %i nsp : %i', t,nspots));
        imgStk=data(:,:,:,1,t);

        % remove spots that don't pass the first round of fitting, anyway
        currentTest = testRatios{t};
        badIdx = currentTest(:,2) < dataProperties.amplitudeCutoff;
        cordList(badIdx,:) = [];
        
        while (~isempty(cordList))
            %             discerned = 0;
            %             if length(cordList) >= maxNumCoords
            %                 error(sprintf('too many coords (%i) - Matlab will run out of memory!',length(cordList)));
            %             end
            [spotsidx mask] = discernspots(cordList,size(imgStk),dataProperties);
            %             discerned = 1;
            %             if length(spotsidx) >= maxNumSpots
            %                 error(sprintf('too many spots (%i) - Matlab will run out of memory!',length(spotsidx)));
            %             end
            idxList=find(mask);
            mskData=imgStk(idxList);

            % do the mixture-model fitting
            [numDist,ncordList,ampList,bg,statistics,debugDataFitTest]=...
                fitTest(mskData,cordList(spotsidx,:),idxList,...
                size(imgStk),dataProperties);

            %             % debug: Collect testValues of intensity-ttest
            %             if ~isempty(debugDataFitTest) && isfield(debugDataFitTest,'testValue')
            %                 tv = cat(2,debugDataFitTest.testValue);
            %                 if ~isfield(debugData,'testValue')
            %                     debugData.testValue = cell(tsteps,1);
            %                 end
            %                 debugData.testValue{t} = ...
            %                     [debugData.testValue{t};ones(length(tv),1)*t,tv(:)];
            %             end

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
                patch=stamp3d(imgStk,[2 2 2],round(ncordList(i,:)));
                %change back to standard coords
                cordOut(t).sp(ct).cord=[ncordList(i,2) ncordList(i,1) ncordList(i,3)];
                cordOut(t).sp(ct).amp=ampList(i);
                cordOut(t).sp(ct).bg=bg;
                cordOut(t).sp(ct).mnint=mean(patch(:));
                % overlapping or single spot
                multi_1=numDist-i+1>size(cordList,1);
                cordOut(t).sp(ct).mult=statistics.multi;
                %                    cordOut.sp(ct).parms=statistics.parms;
                ct=ct+1;
                chi=[chi statistics.chi];
            end; % for i= 1:numDist
        end; %  while (~isempty(cordList))


        % test for zero distance
        cordOut(t).parms=statistics.parms;

        if ~isempty(Q)
            % normally, there should not be a problem here anymore - the
            % tags that were close were fitted jointly and tested already.
            % However, there is the unlikely event that we had single tags
            % that were really badly fitted and thus not significantly
            % different from each other
            cordOut(t) = testdistance(cordOut(t),Q,chi,snr,dataProperties);
            % just add qAmp for the moment; don't check for removed
            % spots
            cordOut(t).statistics.qAmp = qAmp;
        else
            % the if and else here  is only because the nspots>=MAXSPOTS
            cordOut(t).parms=[];
        end;

    end; % if isfield(cord(t),'sp') &  ~isempty(cord(t).sp)
    %     catch
    %         % remember length of spotsIdx
    %         if discerned
    %             maxNumSpots = min(maxNumSpots,length(spotsidx));
    %         else
    %             maxNumCoords = min(maxNumCoords,length(cordList));
    %         end
    %
    %         disp(sprintf('MMF aborted in frame %i:\n%s',t,lasterr))
    %     end
    % progressbar
    if verbose
        mywaitbar(t/tsteps/2,h,tsteps*2);
    end
end; % loop time

warning(warningState)

if verbose
    close(h);
end
