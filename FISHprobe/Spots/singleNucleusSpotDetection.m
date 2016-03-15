function [nucleiStruc, dataProperties] = singleNucleusSpotDetection(metaData, nucleiStruc, varargin)
%singleNucleusSpotDetection detects spots from each multi-channel 3D stack
%of single nucleus
%   Detailed explanation goes here

% To be optimized
% metaData = handles.data.metadata; (Multi-channel 3D stack information)
% handles is a parameter from InvivoCytometer_2.0_source_code/code_package/CellSegmentationQualityAnnotator.m

% 03/2016 Ning

p = inputParser;
p.CaseSensitive = false;
p.addRequired( 'imInput', @(x) ( ismember( ndims(x), [2,3] ) ) );
p.parse( imInput );

% Need to figure out how to read image file and extract property
% information

% imPath = metaData.dataFilePath;

dataProperties.PIXELSIZE_XY = metaData.pixelSize(1);
if metaData.pixelSize(1) ~= metaData.pixelSize(2)
    error('Pixel size on X and Y are not identical.')
end
dataProperties.PIXELSIZE_Z = metaData.pixelSize(3);

% ============================================================
% Define dataProperties parameters
% Generate dataProperties.FILTERPRM
% if isempty(nucleiStruc.numAperture_)
%     dataProperties.NA = input('Enter Numerical Aperture > ');
% else
%     dataProperties.NA=nucleiStruc.numAperture_;
% end
% 
% if isempty(nucleiStruc.channels_.emissionWavelength_)
%     dataProperties.WVL = input('Enter Emission Wavelenth in um > ');
% else
%     dataProperties.WVL = nucleiStruc.channels_.emissionWavelength_/1000;
% end

dataProperties.NA=1.4000;
dataProperties.WVL=0.5250;
dataProperties.refractiveIndex = 1.51;

% Debug mode
% dataProperties.NA = input('Enter Numerical Aperture > ');
% dataProperties.WVL = input('Enter Emission Wavelenth in um > ');
% lenseType = input('Enter lense type (air, water or oil) > ', 's');
% switch lenseType
%     case 'air'
%         dataProperties.refractiveIndex = 1;
%     case 'water'
%         dataProperties.refractiveIndex = 1.33;
%     case 'oil'
%         dataProperties.refractiveIndex = 1.51;
% end

% sigmaCorrection defined by default
dataProperties.sigmaCorrection=[1 1];

% calcFilterParms generates psf size, which is used to define patchsize and
% involved in gaussian filter
% Code borrowed and modified from 
% /home2/nzhang/matlab/applications/FISHprobe/Spots/detect/defaultDataProperties.m

[FT_XY, FT_Z] = calcFilterParms(dataProperties.WVL, ...
    dataProperties.NA, dataProperties.refractiveIndex, 'gauss',...
    dataProperties.sigmaCorrection, ...
    [dataProperties.PIXELSIZE_XY, dataProperties.PIXELSIZE_Z]);

patchXYZ=roundOddOrEven(4*[FT_XY FT_XY FT_Z], 'odd', 'inf');
dataProperties.FILTERPRM = [FT_XY, FT_XY, FT_Z, patchXYZ];
dataProperties.FT_SIGMA = [FT_XY, FT_XY, FT_Z];

Cha = input('Enter the channel for spots detection (red or green) > ', 's');
for nucNum = 1:size(nucleiStruc, 2)
    nuc.dapi = nucleiStruc(nucNum).dapi;
    nuc.red = nucleiStruc(nucNum).red;
    nuc.green = nucleiStruc(nucNum).green;
    
    % 3D image normalization for specific channels and gaussian filter, verbose=0, no waitbar  
    switch Cha
        case 'red'
            nucStack = nuc.red;
        case 'green'
            nucStack = nuc.green;
    end    
    nucStack = filtermovie(nucStack, dataProperties.FILTERPRM, 0);

    spots = spotFindSingleNuc(nucStack, dataProperties);
    nucleiStruc(nucNum).spot = spots.sp;
    
%     spotsPlot3(nucleiStruc, singleChannel3D)
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
        
    % Take MAXNUMSPOTS plus a few spots - we want to be sure that we don't
    % accidentially throw out a good spot, and we need a few bad apples to
    % make the amplitude cutoff work fine. We take between 2 and 10 more
    % spots, depending on MAXNUMSPOTS
    
%     additionalSpots = dataProperties.MAXSPOTS * 0.3;
%     additionalSpots = ceil(additionalSpots);
%     additionalSpots = max(additionalSpots,3);
%     additionalSpots = min(additionalSpots,10);
%     numberOfSpots = dataProperties.MAXSPOTS + additionalSpots;

    % Choose qualified spots number based on mnp value
    [mnpSorted,sortIdx] = sort(mnp(1:ct-1,t),1,'descend');
    
    % Plot cumulative histogram for mnp
    LM = zeros(size(mnpSorted,1),1);
    for i = 1:size(LM)
        LM(i) = size(mnpSorted,1)-i+1;
    end
    plot(mnpSorted,LM,'r*');
    
    % Need to optimize spots selection criteria!!!
    
%     mnpThreshold = input('Enter spottiness threshold > ');
    medianScore = median(mnpSorted);
    distToMedian = mnpSorted-medianScore;
    distThreshold = (mnpSorted(1)-medianScore)/5;
    
    cps = sortIdx(distToMedian-distThreshold>0);
    
    
%     mnpThreshold = 0.001;
%     
%     for qualifiedNum = 1:size(mnpSorted,1)
%         if mnpSorted(qualifiedNum) < mnpThreshold
%             qualifiedNum = qualifiedNum - 1;
%             break
%         else
%             qualifiedNum = qualifiedNum + 1;
%         end
%     end
%     numberOfSpots = qualifiedNum;
%     % cut at either MAXSPOTS+1 or how many we have if it's less
%     cps = sortIdx(1:min(numberOfSpots,length(sortIdx)));
    
    
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
%     spots(t).mnint=mn(t);

        
    % clean memory
    clear FXX FXY FXZ FYX FYY FYZ FZX FZY FZZ
    
end;
