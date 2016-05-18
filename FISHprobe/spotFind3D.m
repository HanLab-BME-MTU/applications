function spots = spotFind3D(fImg, chaParams, detectionMethod, mannualAdjMode)
%spotFindSingleNuc locates fluorescent tags in 3D data
%
% The creteria of spots selection is based on mnp thresholding
% SYNOPSIS cord = spotfind(img)
%
% INPUT img: stack time series
%
% OUTPUT spots: nTimepoints-by-1 structure with fields
%               .sp  structure with fields
%               .cord coordinates
%               .mnint spottiness

% 05/2016 Ning Zhang
% Refer to spotfind.m by dT

patchSize = chaParams.filterParms(4:6);
% init vars
d = floor(patchSize/2);
% inTestD = floor(FILTERSIZE/2); %number of pixels a spot has to be away from the border to be accepted
inTestD = [3,3,1];
movieSize = size(fImg);


%preassign mnp. provide for 100 points
mnpRows = 100;
mnpRowIncrement = 100;
mnp = zeros(mnpRows,1);

% preassign spots
spots = struct('sp',[]);

%intialize counter
ct = 1;
mnplist = [];    %'spottiness'
lst = [];        % list of local maxs
k = [];          % curvature of local maxs

pt = fImg(:,:,:,1);

%norm to 0..100
% Normalization to bitdepth??
pt = 100*pt/max(pt(:));

%find all local max
b = loc_max3Df(fImg(:,:,:,1), [3 3 3]);

[FXX,FXY,FXZ,FYX,FYY,FYZ,FZX,FZY,FZZ] = hessian(pt); % hessian matrix of full intensity dist.

%loop through all local maxs
for i = 1:size(b,1)
    %ignore pixels close to border
    if(all((b(i,:)-inTestD) > 0) && all((b(i,:)+inTestD) <= movieSize(1:3)))
        
        %cut pixels belonging to this local maximum
        patch = stamp3d(pt,patchSize,b(i,:),0);
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
                mnp = zeros(newMnpRows,1);
                mnp(1:mnpRows,:) = mnpTmp;
                mnpRows = newMnpRows;
                
                clear mnpTmp newMnpRows
                
            end
            
            % mnp, the spottiness criterion, is the product of
            % curvature and mean intensity of the local mask.
            % The cutoff might be a bit nicer if we transformed
            % curvature and intensity to [0,1] first, but I leave it
            % for the moment.
            
            switch lower(detectionMethod)
                case 'mnp'
                    mnp(ct) = -k(ct)*mean(patch(:));
                case 'intensity'
                    mnp(ct) = mean(patch(:));
            end
            centp(ct,:) = centroid3D(patch);
            lm(ct,:) = b(i,:);
            ct = ct+1;
            
            
        end;
    end;
end;

% Choose qualified spots number based on mnp value
[mnpSorted,sortIdx] = sort(mnp(1:ct-1),1,'descend');

% Need to optimize spots selection criteria!!!
distThreshold = mnpSorted(1)/5;

% Plot cumulative histogram for mnp
figure
localMax = zeros(size(mnpSorted, 1), 1);
for i = 1:size(localMax)
    localMax(i) = size(mnpSorted, 1)-i+1;
end
plot(mnpSorted, localMax, 'r*');
hold on;
yAxis = ylim;
plot(distThreshold*ones(1,2), [yAxis(1), yAxis(2)], 'b-', 'LineWidth', 2.0)
xlabel('Spottiness Score');
ylabel('Number of Local Maxima');
title('Cumulative Histongram of Spottiness');

if mannualAdjMode
    disp('Please click to choose the spottiness threshold');
    xyThresh = ginput(1);
    distThreshold = xyThresh(1);
    plot(xyThresh(1)*ones(1,2), [yAxis(1), yAxis(2)], 'g-', 'LineWidth', 2.0)
end

hold off;
% close

cps = sortIdx(mnpSorted - distThreshold>0);


if cps  ~=0
    lst = [lm(cps,2) lm(cps,1) lm(cps,3)]-ones(length(cps),1)*(d+1)+centp(cps,:);
    mnplist = mnp(cps);
    j = 1;
    for i=1:size(lst,1)
        % Set boundary for spots center
        if lst(i,:) > patchSize/3
            if lst(i,:) < [size(fImg,2),size(fImg,1),size(fImg,3)]-patchSize/3
                % store coordinates and spottiness
                spots.sp(j).cord = lst(i,:);
                spots.sp(j).mnint = mnplist(i);
                j = j+1;
            end
        end
    end;
end

% clear memory
clear FXX FXY FXZ FYX FYY FYZ FZX FZY FZZ

end