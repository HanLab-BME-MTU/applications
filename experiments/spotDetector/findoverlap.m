function cordOut=findoverlap(data,cord,dataProperties,verbose)
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
if nargin < 3 || isempty(verbose)
    verbose = 1;
end

%CONST DEFINITIONS
MAXSPOTS=dataProperties.MAXSPOTS;

maxNumSpots = 99999; % maximum number of spots that can be fitted simultaneously
maxNumCoords = 99999; % maximum number of coordinates that can be checked by discernspots
tsteps=size(data,5);
if verbose
    h= mywaitbar(0,[],tsteps,'Checking for multiple spots...');
end

% init cordOut. Code changed to no longer overwrite input cord, but to
% build a new structure from scratch. This could potentially be a problem
% with huge structures, but it is much less error-prone
cordOut(1:tsteps) = struct('sp',[],'mnint',[],'statistics',[],'parms',[]);

for t=1:tsteps
    Q=[];
    chi=[];
    snr=[];
    ct=1;
    cord(t).statistics=[];
    cord(t).parms=[];

    %======DEBUG
    %         t;
    %         if any(t == [999])
    %             disp('press dbstep to continue!')
    %             keyboard
    %         end
    %======DEBUG

    % if there are too many local maxima, the code could crash
    % due to memory issues. therefore, try/catch to have only one bad
    % frame, not lose an entire movie

    try
        % if no spots in this time frame -> go back
        if isfield(cord(t),'sp') &  ~isempty(cord(t).sp)
            cordList=cat(1,cord(t).sp.cord);
            %change to matlab coords
            tc2=cordList(:,2);
            cordList(:,2)=cordList(:,1);
            cordList(:,1)=tc2;

            nspots=size(cordList,1);

            imgStk=data(:,:,:,1,t);
            %resize data if necessary (odd size)
            idx=~rem(size(imgStk),2);
            imgStk=imgStk(1:end-idx(1),1:end-idx(2),1:end-idx(3));

            while (~isempty(cordList))
                discerned = 0;
                if length(cordList) >= maxNumCoords
                    error(sprintf('too many coords (%i) - Matlab will run out of memory!',length(cordList)));
                end
                [spotsidx mask] = discernspots(cordList,size(imgStk),dataProperties);
                discerned = 1;
                if length(spotsidx) >= maxNumSpots
                    error(sprintf('too many spots (%i) - Matlab will run out of memory!',length(spotsidx)));
                end
                idxList=find(mask);
                %             % less than 1.5*MAXSPOTS spots found -- do not do that
                %             anymore - fitTest should be able to remove the respective tags!
                %             if nspots<1.5*MAXSPOTS
                mskData=imgStk(idxList);
                [numDist,ncordList,ampList,bg,statistics]=fitTest(mskData,cordList(spotsidx,:),idxList,size(imgStk),dataProperties);
                %             else
                %                 numDist=length(spotsidx);
                %                 ncordList=cordList(spotsidx,:);
                %                 ampList=zeros(numDist,1);
                %                 bg=0;
                if isempty(ncordList)
                    statistics.parms=[];
                    statistics.multi=0;
                    statistics.chi=[];
                    statistics.Qxx=[];
                    statistics.snr=[];
                end;
                %define new cordLst
                tempidx=1:size(cordList,1);
                tempidx(spotsidx)=0;
                tempidx=nonzeros(tempidx);
                cordList=cordList(tempidx,:);

                % add Q of spot(s)
                Q=blkdiag(Q,statistics.Qxx);
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
            cordOut(t).mnint=cord(t).mnint;
            if ~isempty(Q)
                cordOut(t) = testdistance(cordOut(t),Q,chi,snr,dataProperties);
            else
                % the if and else here  is only because the nspots>=MAXSPOTS
                cordOut(t).parms=[];
            end;

        end; % if isfield(cord(t),'sp') &  ~isempty(cord(t).sp)
    catch
        % remember length of spotsIdx
        if discerned
            maxNumSpots = min(maxNumSpots,length(spotsidx));
        else
            maxNumCoords = min(maxNumCoords,length(cordList));
        end

        disp(sprintf('MMF aborted in frame %i:\n%s',t,lasterr))
    end
    % progressbar
    if verbose
        mywaitbar(t/tsteps,h,tsteps);
    end
end; % loop time

if verbose
    close(h);
end
