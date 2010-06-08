function [data] = determinePitDensities(data)
% determinePitDensities calculates density of pits in whole movie
% 
% INPUT:    data    =   struct containing all data for one condition (e.g. 
%                       clathrin control)
%           
%
% OUTPUT:   data    =   structure with added field
%                       .pitDensity 
%                       which is a vector containing the density of objects
%                       in each frame of the movie
%
%           NOTE: the dimension of the density is 
%           # of objects per square pixel
%
%           NOTE2: An earlier version used convex hull to determine the
%           area; this can be tricky for more complicated cell shaped, so
%           this version determines the area using the overlaid and
%           filtered positions of all detected objects
%
%
% Dinah Loerke, last changed 03/19/2008
%               modified 07/29/2008

localDensity = zeros(1, length(data));

for i = 1:length(data)
    
    trackInfoPath = [data(i).source 'TrackInfoMatrices' filesep 'trackInfo.mat'];
    
    if (exist(trackInfoPath, 'file')==2)
        ti = load(trackInfoPath);
        if isfield(ti,'trackInfo')
            currTrackInfo = full(ti.trackInfo);
        elseif isfield(ti,'trackInfoMat')
            currTrackInfo = full(ti.trackInfoMat);
        end
        
        nFrames = size(currTrackInfo, 2)/8;
        
        % extract all x and y positions of object
        cmx = currTrackInfo(:,1:8:end);
        cmy = currTrackInfo(:,2:8:end);
        
        % number of objects detected in each frame
        numPitsFrame = zeros(1,nFrames);
        for k = 1:nFrames
            numPitsFrame(k) = sum(cmx(:,k)>0);
        end
        
        % all positions together
        allx = nonzeros(cmx(:));
        ally = nonzeros(cmy(:));
        %nx = length(allx(:));
        %ny = length(ally(:));
        %numPits = nx/nFrames;
        
        % earlier version used convex hull; this can be tricky for more
        % complicated cell shaped, so allow more complex edge
        % calculate 'basis mask', which is the overlay of all detected objects
        % first, set central pixels of detected objects to 1
        mask_detection = zeros( round(max(allx)), round(max(ally))  );
        detxpos = max(1, round(allx));
        detypos = max(1, round(ally));
        mask_detection(sub2ind(size(mask_detection), detxpos, detypos)) = 1;

        % second, dilate the single pixels to little discs
        mask_detectdilate = imdilate(logical(mask_detection), strel('disk',5));
        
        %figure; imagesc(mask_detectdilate); colormap(gray(256)); axis image;
        
        % third, fill the resulting area and filter with Gaussian profile to
        % create a continuous area
        mask_filled = imfill(double(mask_detectdilate), 'holes');
        mask_close = imclose(mask_filled, strel('disk', 20));
        
        totalarea = sum(mask_close(:));
        denPitsFrame = numPitsFrame/totalarea;        
              
        data(i).pitDensity = denPitsFrame;
        localDensity(i) = nanmean(denPitsFrame);
    else
        data(i).pitDensity = NaN;
    end 
end
fprintf('Pit density %2.4f +- %2.4f obj. per square pixel.\n', nanmean(localDensity), nanstd(localDensity));   