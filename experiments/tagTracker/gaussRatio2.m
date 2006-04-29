function [infoStruct] = gaussRatio(infoStruct,nCoords,constants,movieFrame)
 
% to calculate the ratioMap, we will evaluate the Gauss of every Tag
% at every gridpoint from which we'll later extract intensities.
 
nTags = length(infoStruct);
% catenate coordList
allCoords = cat(1,infoStruct.coordList);
% coords for the first Tag will be allCoords(allCoordsIdx(1):aci(2)-1)
allCoordsIdx = [cumsum(nCoords)];
allGauss = zeros(allCoordsIdx(end)-1, nTags);
 
for iTag = 1:nTags
 
    % get Gauss values for every center at every coordinate position
    allGauss(:,iTag) = GaussListND(allCoords,...
        constants.filterSigma,...
        infoStruct(iTag).centerCoord) * infoStruct(iTag).amp;
 
end

% set minimum intensity of Gauss to 1e-4*min(amp). Check later about possible
% improvements
allGauss = max(allGauss,1e-4*min(cat(1,infoStruct.amp)));
 
% now calculate the ratioMap
for iTag = 1:nTags
 
    % get indexList into allGauss
    idxList = [allCoordsIdx(iTag):allCoordsIdx(iTag+1)-1]';
    
    % read coordList
    coordList = infoStruct(iTag).coordList;
 
    % Gauss one col of allGauss, ratio is this col divided by the sum
    % of all cols (for the idxList)
    infoStruct(iTag).gaussList = allGauss(idxList,iTag);
 
    % currently we divide by the sum of ALL the distributions. In the
    % origninal code, there was no division below a certain threshold
    infoStruct(iTag).ratioList = infoStruct(iTag).gaussList./...
        sum(allGauss(idxList,:),2); % could put sum outside the loop
 
    % with the ratioMap, we can go and read out intensities from the
    % movie frames (which had their background removed in the function
    % call)
    %!!! interpolation uses the image coordinate system!!!
    
    feature('accel','off');  %TMW added
    if ~strcmp(constants.interpolation{1},'*cubic')
        feature('accel','on');   %TMW added
        infoStruct(iTag).intList = interp3(movieFrame,...
            coordList(:,2),...
            coordList(:,1),...
            coordList(:,3),constants.interpolation{1}).*...
            infoStruct(iTag).ratioList;
    else
        % use cubic interpolation directly
        [nrows,ncols,npages] = size(movieFrame);
        [s,t,w] = deal(zeros(size(coordList,1),1));
        
        s(:) = coordList(:,2);
        w(:) = coordList(:,3);
        t(:) = coordList(:,1);
               
        % steps sizes for indexing into a vectorized 3D array
        nw = (nrows+2)*(ncols+2);
        ndx = floor(t)+floor(s-1)*(nrows+2)+floor(w-1)*nw;
 
        % get sub-pixel positions
        feature('accel','off');  %TMW added
        d = find(s==ncols);
        feature('accel','on');   %TMW added
        s(:) = (s - floor(s));
        if length(d)>0, s(d) = s(d)+1; ndx(d) = ndx(d)-nrows-2; end
        d = find(t==nrows);
        t(:) = (t - floor(t));
        if length(d)>0, t(d) = t(d)+1; ndx(d) = ndx(d)-1; end
        d = find(w==npages);
        w(:) = (w - floor(w));
        if length(d)>0, w(d) = w(d)+1; ndx(d) = ndx(d)-nw; end

        % Expand image - use this potentially for gauss-filtering, too!.
            vv = zeros(size(movieFrame)+2);
            vv(2:nrows+1,2:ncols+1,2:npages+1) = movieFrame;
        vv(1,:,:)        = 3*vv(2,:,:)       -3*vv(3,:,:)     +vv(4,:,:); % Y edges
        vv(nrows+2,:,:)  = 3*vv(nrows+1,:,:) -3*vv(nrows,:,:) +vv(nrows-1,:,:);
        vv(:,1,:)        = 3*vv(:,2,:)       -3*vv(:,3,:)     +vv(:,4,:); % X edges
        vv(:,ncols+2,:)  = 3*vv(:,ncols+1,:) -3*vv(:,ncols,:) +vv(:,ncols-1,:);
        vv(:,:,1)        = 3*vv(:,:,2)       -3*vv(:,:,3)     +vv(:,:,4); % Z edges
        vv(:,:,npages+2) = 3*vv(:,:,npages+1)-3*vv(:,:,npages)+vv(:,:,npages-1);
        nrows = nrows+2; ncols = ncols+2; npages = npages+2;
 
        % Interpolate.
        F = zeros(size(s));
        for iw = 0:3,
 
            switch iw
                case 0
                    ww = ((2-w).*w-1).*w;
                case 1
                    ww = (3*w-5).*w.*w+2;
                case 2
                    ww = ((4-3*w).*w+1).*w;
                case 3
                    ww = (w-1).*w.*w;
            end
            for is = 0:3,
                switch is
                    case 0
                        ss = ((2-s).*s-1).*s;
                    case 1
                        ss = (3*s-5).*s.*s+2;
                    case 2
                        ss = ((4-3*s).*s+1).*s;
                    case 3
                        ss = (s-1).*s.*s;
                end
                for it = 0:3,
                    switch it
                        case 0
                            tt = ((2-t).*t-1).*t;
                        case 1
                            tt = (3*t-5).*t.*t+2;
                        case 2
                            tt = ((4-3*t).*t+1).*t;
                        case 3
                            tt = (t-1).*t.*t;
                    end
                    F(:) = F + vv(ndx+(it+is*nrows+iw*nw)).*ss.*tt.*ww;
                end
            end
        end
        F(:) = F/8;
 
% assign output
        infoStruct(iTag).intList = F(:).*...
            infoStruct(iTag).ratioList;
        
        
    end
 
end


