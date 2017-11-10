function [ D, rD ] = analyzeLaminsNPCDistance( MD, lamin_z_position, npc_z_position )
%analyzeLaminsNPCDistance Analyze distance between lamin meshwork

if(nargin < 1 || isempty(MD))
    MD = MovieData.load;
elseif(ischar(MD))
    MD = MovieData.load(MD);
end

cd(MD.movieDataPath_);
mkdir('analyzeLaminsNPCDistance');
cd('analyzeLaminsNPCDistance');

if(nargin < 3)
    if(exist('z_parameters.mat','file'))
        load('z_parameters.mat');
        lamin_z_position
        npc_z_position
    else
        hmv = movieViewer(MD);
        defaultLaminPos = round(MD.zSize_/2);
        defaultNPCPos = defaultLaminPos-2;
        zpos = inputdlg({'Lamin Z Position','NPC Z Position'}, ...
            'Enter Z Positions for Lamins and NPCS', ...
            1, ...
            {num2str(defaultLaminPos),num2str(defaultNPCPos)}, ...
            struct('WindowStyle','normal'));
        lamin_z_position = str2double(zpos{1});
        npc_z_position = str2double(zpos{2});
        close(hmv);
    end
end

save('z_parameters.mat','lamin_z_position','npc_z_position');

gcp

slash_location = strfind(MD.movieDataPath_,'\');
md_title = MD.movieDataPath_(slash_location(end-1):end);

% Lamin image
I = MD.channels_(1).loadImage(1,lamin_z_position);
I_npc = MD.channels_(2).loadImage(1,npc_z_position);

F = OrientationSpaceRidgeFilter(1./2/pi./2,[],8,'none');
R = I*F;
R3 = R.getResponseAtOrderFT(3,2);

import lamins.functions.*;
import lamins.classes.*;

% Get NLMS
original.maxima = R.getRidgeOrientationLocalMaxima;
[nlms,nlms_offset] = nonLocalMaximaSuppressionPrecise(R3.a,original.maxima);

% Get Lamin Threshold
nlms_nz = nlms(nlms ~= 0 & ~isnan(nlms));
% Perhaps detect when there are too many inliers
[~,inliers] = detectOutliers(nlms_nz);
if(length(inliers)./length(nlms_nz) > 0.85)
    T = thresholdOtsu(nlms_nz(inliers));
else
    T = thresholdOtsu(nlms_nz(:));
end
nlms_offset(nlms < T) = NaN;

nlms_mip = nanmax(nlms,[],3);
hfig = figure;
imshow(nlms_mip,[]);
title(md_title);
saveas(hfig,'Lamin_nlms_mip.fig');
saveas(hfig,'Lamin_nlms_mip.png');

% Get Mask
% M = imfill(nlms_mip > T,'holes');
% figure; imshow(M,[]);
% M = imopen(M,strel('disk',5));
% cc = bwconncomp(M);
% [~,idx] = max(cellfun('prodofsize',cc.PixelIdxList));
% mask = labelmatrix(cc) == idx;
% figure; imshowpair(I,mask);
mask = maskFromSteerable(R3);
hfig = figure;
imshowpair(I,mask);
title(md_title);
saveas(hfig,'Nucleus_mask.fig');
saveas(hfig,'Nucleus_mask.png');

% Apply mask to NLMS
nlms_offset(~repmat(mask,[1 1 size(nlms_offset,3)])) = NaN;

[X,Y] = meshgrid(1:size(nlms,2),1:size(nlms,1));

% Get sub-pixel NLMS points
XX = joinColumns(X+cos(original.maxima).*nlms_offset);
YY = joinColumns(Y+sin(original.maxima).*nlms_offset);

% Find sub-pixel NLMS angles
s = ~isnan(XX);
A2 = squeeze(interp3(R.a,repmat(XX(s),1,1,17),repmat(YY(s),1,1,17),ones(size(YY(s))).*shiftdim(1:17,-1)));
[maxima2] = interpft_extrema(A2,2);
[~,minidx] = min(min(abs(maxima2 - original.maxima(s)*2),2*pi-abs(maxima2 - original.maxima(s)*2)),[],2);
minind = sub2ind(size(maxima2),(1:length(minidx)).',minidx);
MM = maxima2(minind)/2;

% Extrapolate additional points from sub-pixels points to decrease
% discretization artifact
LL = [XX(s) YY(s);
      XX(s)+0.5*cos(MM+pi/2) YY(s)+0.5*sin(MM+pi/2); XX(s)-0.5*cos(MM+pi/2) YY(s)-0.5*sin(MM+pi/2);
      XX(s)+0.25*cos(MM+pi/2) YY(s)+0.25*sin(MM+pi/2); XX(s)-0.25*cos(MM+pi/2) YY(s)-0.25*sin(MM+pi/2)];

% Find NPCs
psd = pointSourceDetection(I_npc,2);
[~,psd.inliers] = detectOutliers(psd.A);
if(length(psd.inliers)./length(psd.A) > 0.8)
    psd.T = thresholdRosin(psd.A(psd.inliers));
else
    psd.T = thresholdRosin(psd.A(:));  
end
psd.s = psd.A > psd.T;
psd.s2 = mask(sub2ind(size(I_npc),round(psd.y),round(psd.x)));
psd.s3 = psd.s&psd.s2;

% Plot NPCs
hfig = figure;
imshow(I_npc,[]);
hold on;
scatter(psd.x(psd.s3),psd.y(psd.s3));
title(md_title);
saveas(hfig,'NPC_detection.fig');
saveas(hfig,'NPC_detection.png');

% Find Distance between lamin and NPC
[~,D] = knnsearch(LL,[psd.x(psd.s3).' psd.y(psd.s3).']);

% Find random distance distribution
rxy = rand(60000,2)*1024;
maskch = regionprops(mask,'ConvexHull');
rxyip = inpolygon(rxy(:,1),rxy(:,2),maskch.ConvexHull(:,1),maskch.ConvexHull(:,2));
rxy = rxy(rxyip,:);
[~,rD] = knnsearch(LL,rxy);

% Pair histogram figure
hfig = figure; h = histogram(D.*MD.pixelSize_,0:10:(max(D.*MD.pixelSize_)+10),'Normalization','probability');
hold on; h2 = histogram(rD.*MD.pixelSize_,h.BinEdges,'Normalization','probability');
xlabel('NPC - Lamin Nearest Neighbor Distance (nm)');
ylabel('Probability');
title(md_title);
saveas(hfig,'NPC-Lamin_Distance_Pair_Histogram.fig');
saveas(hfig,'NPC-Lamin_Distance_Pair_Histogram.png');


% Delta Histogram Figure
hfig = figure; bar(h.BinEdges(1:end-1)+5,h.Values - h2.Values);
xlabel('NPC - Lamin Nearest Neighbor Distance (nm)');
ylabel('Probability Difference From Random');
ylim([-0.04 0.04]);
title(md_title);
saveas(hfig,'NPC-Lamin_Distance_DifferenceFromRandom.fig');
saveas(hfig,'NPC-Lamin_Distance_DifferenceFromRandom.png');

% Show Lamin and NPC colocation
hfig = figure;
imshowpair(I,imadjust(I_npc));
hold on; plot(LL(:,1),LL(:,2),'y.');
hold on; scatter(psd.x(psd.s3),psd.y(psd.s3),'c.');
title(md_title);
saveas(hfig,'NPC-Lamin_Detection.fig');
saveas(hfig,'NPC-Lamin_Detection.png');

keyboard;

end

