function [results] = CorrelateData2Pos_intensity(positions, data, rdist);
% CorrelateData2Pos_intensity correlates the intensity in the data (image 
% stack) to the reference positions contained in the mpm file; rdist 
% designates the desired distance around the reference positions
% SYNOPSIS [results] = CorrelateData2Pos_intensity(positions, data, rdist);
% INPUT:    positions:  reference positions, in mpm format
%                       NOTE: time information in the output results (e.g.
%                       intensity in n frames before and after an
%                       internalization event) is encoded in the third
%                       dimension of the mpm file (multiple layers =
%                       multiple time points)
%           data:       data
%                       NOTE: in this function, data needs to be in mpm
%                       format
%           rdist:      relevant distance around reference positions
%                       (can be vector)
%
% OUTPUT:   results:    intensity as a function of distance and time
%                       results has the formats
%                       row 1: rdist
%                       row 2: intensity profile for t=1
%                       row n+1: intensity profile for t=n
%                       where t is the time point specified by the layer of
%                       the positions mpm file
%
%% NOTE1: the function requires that positions is an MPM file; in other 
% versions, zeros in the MPM file are counted as real positions, but this
% isn't possible here because there's no corresponding intensity in the
% image at non-positive index positions; thus zero inputs are converted
% into nan
%% NOTE2 : data needs to be an image stack OR a file containing the 
% subsequent pathnames/filenames of the relevant images; the number of 
% images (or number of specified files) needs to the same as the number of 
% frames in the mpm file
%
% Dinah Loerke, July 24, 2007


%check inputs
[sx,sy,sz] = size(positions);
numf_mpm = round(sy/2);
numt = max(max(positions(:,:,2)));

% determine if data is cell array
if iscell(data)
    datatype = 'cell';
    numf_im = length(data);
else
    datatype = 'stack';
    [sx2,sy2,numf_im] = size(data);
end

numf = min(numf_im,numf_mpm);

% if numf_mpm~=numf_im
%     error('number of frames in positions matrix and number of images don''t match');
% else
%     numf = numf_mpm;
% end

positions(positions==0) = nan;

% initialize results matrix
results = nan*rdist; results(length(results)) = [];

% maximum distance
rad = max(rdist);



% length of the distance vector/dimension
if length(rdist)>1
    rlen = length(rdist)-1;
else
    rlen=1;
end

% initialize results matrix
intensityMat = nan*zeros(sx,numt,rlen);




% (relative) index positions of the pixels that have the specified distance
% from a pixel of position n in an image of this size
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

% loop over all frames
for i=1:numf
    
    fprintf('frame #%03d',i);
    
    %current image
    if strcmp(datatype,'stack')
        currIm = data(:,:,i);
    elseif strcmp(datatype,'cell')
        currIm = imread(data{i});
    else
        error('format not found')
    end
    [ix,iy] = size(currIm);
    
    checkIm = currIm;
             

    % current reference positions
    currx = full(positions(:,2*i-1,1));
    curry = full(positions(:,2*i,1));
    upos = find( isfinite(currx) & isfinite(curry) & (currx>0) & (curry>0) );
    uposx = currx(upos);
    uposy = curry(upos);
    upost = full(positions(upos,2*i,2));


    % loop over all reference points in this frame
    for k=1:length(uposx)

        %fprintf(' obj #%04d',k);


         % initialize
        intensityVector = nan*rdist;


        % pixel positions for this object 
        XobjPos = round(uposx(k));
        YobjPos = round(uposy(k));
        XusePos = XobjPos + XinPix;
        YusePos = YobjPos + YinPix;
        % NOTE: maybe correct distances for integer-value rounding
        % effects????
        DusePos = DinPix;

        % if the object is close to the image edge, check that none of the 
        % pixels are outside the image; note x-y switch for image dim
        if ( (min(XobjPos,abs(iy-XobjPos))<=rad) | (min(YobjPos,abs(ix-YobjPos))<=rad) )
            goodPos = find( (XusePos>0) & (YusePos>0) & (XusePos<iy) & (YusePos<ix) );
            XusePos = XusePos(goodPos);
            YusePos = YusePos(goodPos);
            DusePos = DinPix(goodPos);
        end

                        
        % these are the x-y positions for the reference coordinate
        % system; to match up with the image coordinate system, x and y
        % have to be switched
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
        intensityMat(upos(k),upost(k),:) = intensityVector;


    end % of loop over k (individual objects in frame)

     
    %end % of loop over t (layers of mpm)
    fprintf('\b\b\b\b\b\b\b\b\b\b');
    
    
end % of loop over i (frames)
fprintf('\n');

% % result is the intensity profile average for all frames in the movie
% resultsTotal = nanmean(intensityProfile_movie,1);
% 
% % shuffle the third dimension into the first
% for n=1:sz
%     results(n,:) = resultsTotal(1,:,n);
% end
% % add used rdist vector to the first row for reference
% results = [rdist(1:(length(rdist)-1)); results];

results = intensityMat;

end % of function





            
    