function [beadDistAll,beadDensity] = estimateBeadDistance(img,pixSize,sigma,h)
%function [beadDistAll,beadDensity] = estimateBeadDistance(img) will identify
%beads, and quantify statistics about bead-to-bead distance and bead density.
%   input:  
%           img:                bead image
%           pixSize:            pixel size in nm (output from MovieData)
%           sigma:              point source sigma (from MD)
%           h:                  header for the figure
%   output:
%           beadDistAll:        cell array containing closest distance (in
%                               micron) of all individual beads
%           beadDensity:        the number of beads per area (#/um2)
% May 2020, Sangyoon Han

% sigma =getGaussianPSFsigmaFromData(img);
if nargin<4
    h=[];
end
%% setting
pixSizeUm = pixSize/1000; %0.090; %um
maxDistUm = 1.5;
maxDist= maxDistUm/pixSizeUm;
%% detection
pstruct = pointSourceDetection(img,sigma,'FitMixtures',true,'MaxMixtures',3);
% hold off
% figure(1), imshow(img,[]), hold on
% plot(pstruct.x,pstruct.y,'ro')
[idx,dist] = KDTreeBallQuery([pstruct.x' pstruct.y'],[pstruct.x' pstruct.y'],maxDist);
numDist = cellfun(@length,dist);
dist2 = dist(numDist>1);
distAll=cellfun(@(x) x(2), dist2);
distAllUm = distAll*pixSizeUm;
%% Drawing bead-to-bead spacing
% idxAll = cellfun(@(x) x(2), idx(numDist>1));
% validPoints = find(numDist>1);
% % plot([pstruct.x(validPoints)' pstruct.x(idxAll(validPoints))'], ...
% %     [pstruct.y(validPoints)' pstruct.y(idxAll(validPoints))'],'r')
% h=figure(1);
% for ii= find(numDist>1)'
%     curPoint = [pstruct.x(ii) pstruct.y(ii)];
%     correspondingPoint = [pstruct.x(idx{ii}(2)) pstruct.y(idx{ii}(2))];
%     plot([curPoint(1) correspondingPoint(1)], [curPoint(2), correspondingPoint(2)],'r')
% end
%% drawing detection
if ~isempty(h)
    imshow(img,[])
    hold on
    plot(pstruct.x, pstruct.y, 'ro')
end
%% Bead density
% There is some boundary where beads are not detected. Excluding those area
curAreaPix = (size(img,1)-ceil(4*sigma)*2) * (size(img,2)-ceil(4*sigma)*2);
curArea = curAreaPix * pixSizeUm^2;
%% Output
beadDistAll = distAllUm; %in um
beadDensity = numel(pstruct.x)/curArea; % in #/um2
% beadDensity = 1./sqrt(mean(beadDistAll));
disp(['Mean distance= ' num2str(mean(distAllUm),3) ' um, Bead density= ' num2str(beadDensity,3) ' /um^2.'])
end

