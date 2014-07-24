function depths = edgeDepths(edgePaths,mask,distX)
%EDGEDEPTHS calculates the depth of each point on each edge within the input mask
% 
% depths = edgeDepths(edgePaths,mask)
% 
% ADD MORE HELP YOU LAZY ASS!!
% 
% Hunter Elliott
% 5/2011
%

if nargin < 3 || isempty(distX)    
    distX = bwdist(~mask);
end 
    
nEdges = numel(edgePaths);

depths = cell(nEdges,1);

for j = 1:nEdges
    
    nPts = size(edgePaths{j},1);
    %TEMP - Replace rounding with linear interpolation??? -HLE
    depths{j} = arrayfun(@(x)(distX(round(edgePaths{j}(x,1)),...
                                       round(edgePaths{j}(x,2)),...
                                       round(edgePaths{j}(x,3)))),...
                                                 1:nPts)';

end