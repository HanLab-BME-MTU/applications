function [iMat] = extractIntensity_MPMfromImageStack(MPM, imageStack, rdist);
% reads out the intensities in the image stack at the reference positions 
% specified in the mpm file; rdist designates the desired distance around 
% the reference positions
%
% SYNOPSIS [iMat] = extractIntensity_MPMfromImageStack(MPM, imageStack,rdist);
% 
% INPUT:    MPM:        mpm file of object positions, in the format
%                       [x1 y1 x2 y2 x3 y3.... xn yn];
%           imageStack: image Stack from which intensity is read
%           rdist:      relevant distance around reference positions
%                       (can be vector)
%
% OUTPUT:   iMat:    
%
% NOTE1: the function requires that positions are in MPM format; in other 
% versions, zeros in the MPM file are counted as real positions, but this
% isn't possible here because there's no corresponding intensity in the
% image at non-positive index positions; thus, zero inputs are converted
% into nan
% NOTE2 : image/parameter data needs to be an image stack OR a file 
% containing the subsequent pathnames/filenames of the relevant images; the 
% number of images (or number of specified files) needs to the same as the 
% number of frames in the mpm file
%
% Dinah Loerke, January 30, 2009



%% =======================================================================
%  check inputs and set defaults
% ========================================================================

[sx,sy,sz] = size(MPM);
numf_mpm = round(sy/2);

% determine if imageStack is cell array
if iscell(imageStack)
    datatype = 'cell';
    numf_im = length(imageStack);
else
    datatype = 'stack';
    [sx2,sy2,numf_im] = size(imageStack);
end

numf = min(numf_im,numf_mpm);

MPM(MPM==0) = nan;

% maximum distance
rad = max(rdist);

% length of the distance vector/dimension
if length(rdist)>1
    rlen = length(rdist)-1;
else
    rlen=1;
end

% initialize results matrix
iMat = nan*zeros(sx,numf,rlen);




%% =======================================================================
%  prepare intensity readout
% ========================================================================

% to read out the intensities within a given distance of the object
% positions, we here pre-define (relative) index positions of the pixels 
% that have the specified distance from a pixel of position n in an image 
% of this size - this allows us to avoid having to run a distance check for
% each individual object

% mini-grid
[miniImX, miniImY] = ndgrid(-ceil(rad):ceil(rad),-ceil(rad):ceil(rad));
miniDist = sqrt(miniImX.^2 + miniImY.^2);
%inside pixels
[XinPix, YinPix] = find(miniDist<=rad);
% subtract ceil(rad)+1 to express pixel positions as relative to a given
% position
XinPix = XinPix-(ceil(rad)+1);
YinPix = YinPix-(ceil(rad)+1);
% corresponding distances of all these points
DinPix = sqrt(XinPix.^2 + YinPix.^2);



%% =======================================================================
%  read intensities
% ========================================================================

% loop over all frames
for i=1:numf
    
    fprintf('frame #%03d',i);
    
    %current image
    if strcmp(datatype,'stack')
        currIm = imageStack(:,:,i);
    elseif strcmp(datatype,'cell')
        %currIm = imread(imageStack{i});
        currIm = importdata(imageStack{i});
    else
        error('format not found')
    end
    [ix,iy] = size(currIm);
    
    checkIm = currIm;
             

    % current valid reference positions
    currx = full(MPM(:,2*i-1,1));
    curry = full(MPM(:,2*i,1));
    upos = find( isfinite(currx) & isfinite(curry) & (currx>0) & (curry>0) );
    uposx = currx(upos);
    uposy = curry(upos);
   
    
    % loop over all reference points in this frame
    for k=1:length(uposx)

        % initialize
        intensityVector = nan*rdist;

        % pixel positions for this object 
        XobjPos = round(uposx(k));
        YobjPos = round(uposy(k));
        XusePos = XobjPos + XinPix;
        YusePos = YobjPos + YinPix;
        % NOTE: maybe in the future correct distances for integer-value 
        % rounding effects????
        DusePos = DinPix;

        % if the object is close to the image edge, check that none of the 
        % pixels are outside the image; note x-y switch for image dim
        if ( (min(XobjPos,abs(iy-XobjPos))<=rad) | (min(YobjPos,abs(ix-YobjPos))<=rad) )
            goodPos = find( (XusePos>0) & (YusePos>0) & (XusePos<iy) & (YusePos<ix) );
            XusePos = XusePos(goodPos);
            YusePos = YusePos(goodPos);
            DusePos = DinPix(goodPos);
        end

        % these are now the x-y positions for the reference coordinate
        % system; to match up with the image coordinate system, x and y
        % have to be switched back
        
        currObjIntensities = [];
        currObjDistances = [];

        for in=1:length(XusePos)
            % write intensities of pixels into matrix
            currObjIntensities(in) = currIm(YusePos(in),XusePos(in));
            % write corresponding distance into matrix
            currObjDistances(in) = DusePos(in);
            % internal control
            % checkIm(YusePos(in),XusePos(in)) = max(checkIm(:));
        end
            
            
            %============display results, uncomment if desired
    %         hold off    
    %         imshow(checkIm,[]);
    %         hold on;
    %         plot(currx,curry,'r.');
    %         pause(0.01);
            %=================================================


        % average intensities for all distance bins
        if length(rdist)>1
            intensityVector = thresholdVector(currObjDistances,currObjIntensities,rdist);
        else
            intensityVector = nanmean(currObjIntensities);
        end

        % the intensity for this objects is entered into the appropriate
        % row of the results matrix, and the columns corresponding to the
        % time shift for this object in this frame
        % NOTE that this format for the results matrix works only because
        % there is only one object/trajectory per row!
        iMat(upos(k),i,:) = intensityVector;


    end % of loop over k (individual objects in frame)

     
    %end % of loop over t (layers of mpm)
    fprintf('\b\b\b\b\b\b\b\b\b\b');
    
    
end % of loop over i (frames)
fprintf('\n');

end % of function





            
    