function spots = spotfind_mammalian(filteredMovieName,dataProperties,verbose,movieLoader)
%SPOTFIND locates fluorescent tags in 3D data
%
% SYNOPSIS cord = spotfind(img)
%
% INPUT img   : stack time series
%       dataProperties: structure with movie properties
%       verbose : (optional) If 1 (default), waitbar is displayed
%       loadedFrames: list of loaded frames for improved waitbar
%
% OUTPUT spots : nTimepoints-by-1 structure with fields
%                .sp  structure with fields
%                   .cord coordinates
%                   .mnint spottiness
%                .COM center of mass of image
%
% This code directly sets a cutoff based on spot intensity

% c: 5/3/01	dT

% optional input arguments
if nargin < 3 || isempty(verbose)
    verbose = 1;
end

if nargin < 4 || isempty(movieLoader)
    if isnumeric(filteredMovieName)
        movieLoader = 'none';
    else
        movieLoader = 'any';
    end
end

%CONST DEFINITIONS
%global PATCHSIZE;
%PATCHSIZE=dataProperties.PATCHSIZE;
FILTERSIZE = dataProperties.FILTERPRM(4:6);
PATCHSIZE = FILTERSIZE;
DEBUG = 0;

% init vars
d=floor(PATCHSIZE/2);
% inTestD = floor(FILTERSIZE/2); %number of pixels a spot has to be away from the border to be accepted
inTestD = [3,3,1];

if strcmp(movieLoader,'none')
    fImg = filteredMovieName;
    loadStruct.loadedFrames = 1:size(fImg,5);
    loadStruct.frames2load = [];

else
    % load filtered movie
    [fImg, movieHeader, loadStruct] = ...
        cdLoadMovie({filteredMovieName,'filtered'}, [], dataProperties);

end

plotTime = round(linspace(1,dataProperties.movieSize(4),20));
fh=figure('Name',dataProperties.name);

% preassign output
% preassign spots
spots(1:dataProperties.movieSize(4),1) = struct('sp',[],'COM',[],...
    'cutoffIntensity',[],'intensityList',[],'coordinateList',[]);
% loop through the movie to collect local maxima
done = 0;
while ~done

    movieSize = size(fImg);
    tsteps=size(fImg,5); % just in case there's only one frame loaded

    if verbose
        if isempty(loadStruct.loadedFrames)
            h= mywaitbar(0,[],tsteps,'Finding spots...');
        else
            h= mywaitbar(0,[],tsteps,...
                sprintf('Finding spots in frames (%i:%i)',...
                loadStruct.loadedFrames(1),loadStruct.loadedFrames(end)));
        end
    end



    %loop through all the loaded points
    for t=1:tsteps
        %preassign mnp. provide for 1000 points.
        % cols: 1: medianIntensity 2-4: centroid within patch 5:7 centerPos of
        % patch
        mnpRows = 1000;
        mnpRowIncrement = 1000;
        mnp2 = zeros(mnpRows,7);
        %intialize counter
        ct=1;
        % current time point:
        mnplist=[];    %'spottiness'
        lst=[];            % list of local maxs
        k=[];             % curvature of local maxs

        pt=fImg(:,:,:,1,t);

        % do not norm, b/c we want to keep the noise level constant
        %norm to 0..100
        %pt=100*pt/max(pt(:));



        %find all local max
        b=loc_max3Df(fImg(:,:,:,1,t),[3 3 3]);

        %[FXX,FXY,FXZ,FYX,FYY,FYZ,FZX,FZY,FZZ]=hessian(pt); % hessian matrix of full intensity dist.

        %loop through all local maxs
        for i=1:size(b,1)
            %ignore pixels close to border
            if(all((b(i,:)-inTestD)>0) & all((b(i,:)+inTestD)...
                    <=[movieSize(1:3)]))

                %cut pixels belonging to this local maximum
                patch = stamp3d(pt,PATCHSIZE,b(i,:),0);
                %             patch=pt(b(i,1)-d(1):b(i,1)+d(1),b(i,2)-d(2):b(i,2)+d(2),b(i,3)-d(3):b(i,3)+d(3));


                %curvature filter
                %k(ct)=curvature3D(patch,[d d d]+1);
                %             k(ct)=det([FXX(b(i,1),b(i,2),b(i,3)) FXY(b(i,1),b(i,2),b(i,3)) FXZ(b(i,1),b(i,2),b(i,3));...
                %                 FYX(b(i,1),b(i,2),b(i,3)) FYY(b(i,1),b(i,2),b(i,3)) FYZ(b(i,1),b(i,2),b(i,3));...
                %                 FZX(b(i,1),b(i,2),b(i,3)) FZY(b(i,1),b(i,2),b(i,3)) FZZ(b(i,1),b(i,2),b(i,3))]);

                % only convex shapes allowed
                %             if k(ct) < 0
                if ct > mnpRows
                    % reassign mnp-matrix
                    %mnpTmp = mnp;
                    newMnpRows = mnpRows + mnpRowIncrement;
                    %                     mnp = zeros(newMnpRows,tsteps);
                    %                     mnp(1:mnpRows,:) = mnpTmp;
                    %mnpRows = newMnpRows;

                    %clear mnpTmp newMnpRows

                    mnpTmp = mnp2;
                    %newMnpRows = mnpRows + mnpRowIncrement;
                    mnp2 = zeros(newMnpRows,7);
                    mnp2(1:mnpRows,:) = mnpTmp;
                    mnpRows = newMnpRows;

                    clear mnpTmp newMnpRows
                end

                % mnp, the spottiness criterion, is the product of
                % curvature and mean intensity of the local mask.
                % The cutoff might be a bit nicer if we transformed
                % curvature and intensity to [0,1] first, but I leave it
                % for the moment.
                %patch2 = stamp3d(pt2,PATCHSIZE,b(i,:),0);
                mnp2(ct,1) = median(patch(:));
                %                 patch2=patch2-min(patch2(:));
                %                 patch2 = patch2/max(patch2(:));
                %                 mnp2(ct,t,2) = std(patch(:));
                %                 mnp2(ct,t,2) = mean(patch(patch(:)<mnp2(ct,t,1)));
                %mnp2(ct,t,2) = median((patch(:)-mnp2(ct,t,1)).^2/numel(patch));
                %                 mnp(ct,t)=-k(ct)*mean(patch(:));
                mnp2(ct,2:4)=centroid3D(patch);
                mnp2(ct,5:7)=b(i,:); % was lm
                ct=ct+1;


            end;
        end;
        
        % remove superfluous rows from mnp2
        mnp2(ct:end,:) = [];

        %figure,plot(squeeze(mnp2(1:ct-1,t,1)),squeeze(mnp2(1:ct-1,t,2)),'.')
        %     figure,subplot(2,2,1),histogram(mnp2(1:ct-1,t,1),'smooth')
        %     subplot(2,2,2),histogram(prod(mnp2(1:ct-1,t,:),3),'smooth');
        %     subplot(2,2,3),histogram(mnp(1:ct-1,t),'smooth');
        %     subplot(2,2,4),plot(squeeze(mnp2(1:ct-1,t,1)),mnp(1:ct-1,t),'.');
        currentData = mnp2(:,1);
        if any(loadStruct.loadedFrames(t)==plotTime)
            figure(fh);
            ah = subplot(4,5,find(loadStruct.loadedFrames(t)==plotTime));
        else
            ah = 0;
        end
%                 figure,%ah=subplot(1,2,1);
%                 ah = gca;
%ah = 0;
        cutValue = cutFirstHistMode_robust(currentData,[0,2],ah);
        %ah=subplot(1,2,2);
        firstCut = currentData>cutValue;
        %         cutoffIntensity = cutFirstHistMode_robust(currentData(firstCut),[0,1],ah);
        cutoffIntensity = cutValue;

        %remember positins etc for the spots that have passed the first
        %cut only
        % -- remember all
        spots(loadStruct.loadedFrames(t)).intensityList = currentData(:,1);
        spots(loadStruct.loadedFrames(t)).coordinateList = ...
            [mnp2(:,6) mnp2(:,5) mnp2(:,7)]-ones(ct-1,1)*(d+1)+mnp2(:,2:4);

        spots(loadStruct.loadedFrames(t)).cutoffIntensity = cutoffIntensity;
        % Take MAXNUMSPOTS plus a few spots - we want to be sure that we don't
        % accidentially throw out a good spot, and we need a few bad apples to
        % make the amplitude cutoff work fine. We take between 2 and 10 more
        % spots, depending on MAXNUMSPOTS

        %     additionalSpots = dataProperties.MAXSPOTS * 0.3;
        %     additionalSpots = ceil(additionalSpots);
        %     additionalSpots = max(additionalSpots,3);
        %     additionalSpots = min(additionalSpots,10);
        %     numberOfSpots = dataProperties.MAXSPOTS + additionalSpots;
        %
        %     [mnpSorted,sortIdx] = sort(mnp(1:ct-1,t),1,'descend');
        %     % cut at either MAXSPOTS+1 or how many we have if it's less
        %     cps = sortIdx(1:min(numberOfSpots,length(sortIdx)));
        %
        %
        %     if cps~=0
        %         lst=[lm(cps,2) lm(cps,1) lm(cps,3)]-ones(length(cps),1)*(d+1)+centp(cps,:);
        %         mnplist=mnp(cps,t);
        %         for i=1:size(lst,1)
        %             % store coordinates and spottiness
        %             spots(t).sp(i).cord=lst(i,:);
        %             spots(t).sp(i).mnint=mnplist(i);
        %         end;
        %     end
        %     spots(t).mnint=mn(t);

        %== ADDED FOR NEW LINKER ==
        % Find the "center of gravity" of the image
        % use int^10 to get good results (maybe we need more for mammalian
        % cells with their bigger frames?)
        spots(loadStruct.loadedFrames(t)).COM = centroid3D(pt,10);


        if verbose
            mywaitbar(t/tsteps,h,tsteps);
        end

        % clean memory
        %clear FXX FXY FXZ FYX FYY FYZ FZX FZY FZZ

    end % loop through loaded timesteps

    % load more
    if isempty(loadStruct.frames2load)
        done = 1;
        if verbose
            close(h);
        end
    else
        % free memory
clear fImg
        [fImg, movieHeader, loadStruct] = ...
            cdLoadMovie(loadStruct.movieType, [], loadStruct);
if verbose
            close(h);
        end
    end % load more

end % while ~done



% plot cutoffIntensity
% ci = cat(1,spots.cutoffIntensity);
% figure,plot(ci);
% % for the overall cutoff: take the average to avoid problems with bimodal
% % distributions
% mci = mean(ci);
% hold on
% plot([1,length(ci)],[mci,mci],'r');

% while the noise intensity remains similar, the signal intensities
% decrease over time. Thus, the cutoff will have to be set lower to catch
% the good spots, and should be calculated via a robust exponential fit.
ci = cat(1,spots.cutoffIntensity);
xFit = robustExponentialFit2(ci,[],verbose);

nSpots = zeros(length(spots),1);

for t = 1:length(spots)
    goodIdx = spots(t).intensityList > xFit(1)*exp(xFit(2)*t);
    
    % fill sp
    spots(t).sp = struct(...
        'cord',mat2cell(spots(t).coordinateList(goodIdx,:),double(goodIdx(goodIdx)),3),...
        'mnint',mat2cell(spots(t).intensityList(goodIdx,:),double(goodIdx(goodIdx)),1));
    nSpots(t) = sum(goodIdx);
    
end

% remove intensity list, coordinateList from spots
spots = rmfield(spots,{'coordinateList','intensityList'});

% plot overall histogram
allInt = catStruct(1,'spots.sp.mnint');
figure('Name',dataProperties.name)
ah=subplot(2,1,1);
histogram(ah,allInt);
subplot(2,1,2);
plot(1:length(spots),nSpots,'.');
        



  
