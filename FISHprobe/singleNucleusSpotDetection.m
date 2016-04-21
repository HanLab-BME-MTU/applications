function [nucleiStruc, dataProperties] = singleNucleusSpotDetection(nucleiStruc, dataProperties, imageData, varargin)
%singleNucleusSpotDetection detects spots from each multi-channel 3D stack
%of single nucleus
%   Detailed explanation goes here

% 04/2016 Ning

% p = inputParser;
% p.CaseSensitive = false;
% p.addRequired( 'imInput', @(x) ( ismember( ndims(x), [2,3] ) ) );
% p.parse();
% parse the input dataProperties!!! (No empty input)


% Define dataProperties parameters
% Generate dataProperties.FILTERPRM
% calcFilterParms generates psf size, which is used to define patchsize and
% involved in gaussian filter
% Code borrowed and modified from 
% /home2/nzhang/matlab/applications/FISHprobe/Spots/detect/defaultDataProperties.m

for chaNum = 1:numel(dataProperties.channel)
    chaName = dataProperties.channel(chaNum).name;
    switch chaName
        case 'dapi'
            continue;
            
        case 'green' 
            greenWVL = dataProperties.channel(chaNum).emissionWavelength;
            [FT_XY, FT_Z] = calcFilterParms(greenWVL, dataProperties.NA, ...
                            dataProperties.refractiveIndex, 'gauss', ...
                            dataProperties.sigmaCorrection, ...
                            [dataProperties.PIXELSIZE_XY, dataProperties.PIXELSIZE_Z]);
            % FT_XY vs psfsigma value from bioformat reader??
            patchXYZ = roundOddOrEven(4*[FT_XY FT_XY FT_Z], 'odd', 'inf');
            dataProperties.channel(chaNum).FILTERPRM = [FT_XY, FT_XY, FT_Z, patchXYZ];
            dataProperties.channel(chaNum).FT_SIGMA = [FT_XY, FT_XY, FT_Z];
            
            for nucNum = 1:numel(nucleiStruc)
                nucStack = nucleiStruc(nucNum).green;
                nucStack = filtermovie(nucStack, dataProperties.channel(chaNum).FILTERPRM, 0);
                
                spots = spotFindSingleNuc(nucStack, dataProperties.channel(chaNum));
                nucleiStruc(nucNum).greenSpot = spots.sp;
                
            end
            
        case 'red'
            redWVL = dataProperties.channel(chaNum).emissionWavelength;
            [FT_XY, FT_Z] = calcFilterParms(redWVL, dataProperties.NA, ...
                dataProperties.refractiveIndex, 'gauss', ...
                dataProperties.sigmaCorrection, ...
                [dataProperties.PIXELSIZE_XY, dataProperties.PIXELSIZE_Z]);
            
            patchXYZ=roundOddOrEven(4*[FT_XY FT_XY FT_Z], 'odd', 'inf');
            dataProperties.channel(chaNum).FILTERPRM = [FT_XY, FT_XY, FT_Z, patchXYZ];
            dataProperties.channel(chaNum).FT_SIGMA = [FT_XY, FT_XY, FT_Z];
            
            for nucNum = 1:numel(nucleiStruc)
                nucStack = nucleiStruc(nucNum).red;
                nucStack = filtermovie(nucStack, dataProperties.channel(chaNum).FILTERPRM, 0);
                
                spots = spotFindSingleNuc(nucStack, dataProperties.channel(chaNum));
                nucleiStruc(nucNum).redSpot = spots.sp;
                
            end
            
        otherwise
            error('Unknown channels detected')
    end
   
end

    spotsPlot3(nucleiStruc, imageData, dataProperties);
    
end


function spots = spotFindSingleNuc(fImg, dataProperties)
%spotFindSingleNuc locates fluorescent tags in 3D data
%
% The creteria of spots selection is based on mnp thresholding
% SYNOPSIS cord = spotfind(img)
%
% INPUT img   : stack time series
%
% OUTPUT spots : nTimepoints-by-1 structure with fields
%                .sp  structure with fields
%                   .cord coordinates
%                   .mnint spottiness

% 03/2016 Ning


FILTERSIZE = dataProperties.FILTERPRM(4:6);
PATCHSIZE = FILTERSIZE;

% init vars
d = floor(PATCHSIZE/2);
% inTestD = floor(FILTERSIZE/2); %number of pixels a spot has to be away from the border to be accepted
inTestD = [3,3,1];
movieSize = size(fImg);
tsteps = size(fImg,5); % just in case there's only one frame loaded


%preassign mnp. provide for 100 points
mnpRows = 100;
mnpRowIncrement = 100;
mnp = zeros(mnpRows,tsteps);

% preassign spots
spots(1:tsteps,1) = struct('sp',[]);

%loop through all time points
for t = 1:tsteps

    %intialize counter
    ct = 1;
    % current time point:
    mnplist = [];    %'spottiness'
    lst = [];        % list of local maxs
    k = [];          % curvature of local maxs

    pt = fImg(:,:,:,1,t);

    %norm to 0..100
    pt = 100*pt/max(pt(:));
    
    %find all local max
    b = loc_max3Df(fImg(:,:,:,1,t),[3 3 3]);

    [FXX,FXY,FXZ,FYX,FYY,FYZ,FZX,FZY,FZZ] = hessian(pt); % hessian matrix of full intensity dist.

    %loop through all local maxs
    for i = 1:size(b,1)
        %ignore pixels close to border
        if(all((b(i,:)-inTestD) > 0) && all((b(i,:)+inTestD) <= movieSize(1:3)))

            %cut pixels belonging to this local maximum
            patch = stamp3d(pt,PATCHSIZE,b(i,:),0);
            % patch=pt(b(i,1)-d(1):b(i,1)+d(1),b(i,2)-d(2):b(i,2)+d(2),b(i,3)-d(3):b(i,3)+d(3));


            %curvature filter
            %k(ct)=curvature3D(patch,[d d d]+1);
            k(ct) = det([FXX(b(i,1),b(i,2),b(i,3)) FXY(b(i,1),b(i,2),b(i,3)) FXZ(b(i,1),b(i,2),b(i,3));...
                FYX(b(i,1),b(i,2),b(i,3)) FYY(b(i,1),b(i,2),b(i,3)) FYZ(b(i,1),b(i,2),b(i,3));...
                FZX(b(i,1),b(i,2),b(i,3)) FZY(b(i,1),b(i,2),b(i,3)) FZZ(b(i,1),b(i,2),b(i,3))]);

            % only convex shapes allowed
            if k(ct) < 0
                if ct > mnpRows
                    % reassign mnp-matrix
                    mnpTmp = mnp;
                    newMnpRows = mnpRows + mnpRowIncrement;
                    mnp = zeros(newMnpRows,tsteps);
                    mnp(1:mnpRows,:) = mnpTmp;
                    mnpRows = newMnpRows;
                    
                    clear mnpTmp newMnpRows
                    
                end
                
                % mnp, the spottiness criterion, is the product of
                % curvature and mean intensity of the local mask.
                % The cutoff might be a bit nicer if we transformed
                % curvature and intensity to [0,1] first, but I leave it
                % for the moment.
                mnp(ct,t) = -k(ct)*mean(patch(:));
                centp(ct,:) = centroid3D(patch);
                lm(ct,:) = b(i,:);
                ct = ct+1;
                
                
            end;
        end;
    end;
        
    % Choose qualified spots number based on mnp value
    [mnpSorted,sortIdx] = sort(mnp(1:ct-1,t),1,'descend');
    
%     % Plot cumulative histogram for mnp
%     LM = zeros(size(mnpSorted,1),1);
%     for i = 1:size(LM)
%         LM(i) = size(mnpSorted,1)-i+1;
%     end
%     plot(mnpSorted,LM,'r*');
    
    % Need to optimize spots selection criteria!!!
    
%     mnpThreshold = input('Enter spottiness threshold > ');
    medianScore = median(mnpSorted);
    distToMedian = mnpSorted-medianScore;
    distThreshold = (mnpSorted(1)-medianScore)/5;
    
    cps = sortIdx(distToMedian-distThreshold>0);
    
    
    if cps  ~=0
        lst = [lm(cps,2) lm(cps,1) lm(cps,3)]-ones(length(cps),1)*(d+1)+centp(cps,:);
        mnplist = mnp(cps,t);
        j = 1;
        for i=1:size(lst,1)
            % Set boundary for spots center
            if lst(i,:) > PATCHSIZE/3
                if lst(i,:) < [size(fImg,2),size(fImg,1),size(fImg,3)]-PATCHSIZE/3
                    % store coordinates and spottiness
                    spots(t).sp(j).cord = lst(i,:);
                    spots(t).sp(j).mnint = mnplist(i);
                    j = j+1;
                end
            end
        end;
    end
        
    % clear memory
    clear FXX FXY FXZ FYX FYY FYZ FZX FZY FZZ
    
end;
end