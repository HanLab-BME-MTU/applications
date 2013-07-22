function [Ouput]=meanShift_Merger_poly(ShiftArr, PointList,OverLap)
% meanShift_Merger_poly takes a cell array with several separate meanShift
%  and merges clusters if their convex huls overlap by Overlap. This is
%  interated until the number of resulting clusters are stable.
%
% Inputs: 
%           ShiftArr, a cell array containing the results of
%                     MeanShiftClustering
%
%           PointList, a cell array containing the pnts used to gerate
%                      ShiftArr. Use the output of exchangePaintAlignment
%           
%           OverLap, the percentage of Overlap required for megering
%                    cluster
%
% Outputs: 
%           Output, A struct containing the merged clusters
%                   
%                   .com  : center of mass of the new cluster
%
%                   .pnts : coordinates of point in a cluster and the
%                           indenty of which data set it came from
%                           [x,y,dx,dy,Id]
%
%
% Writen 7/19/2013 by Jeffrey Werbin, Harvard Medical School
%


n = numel(ShiftArr);

if n ~= numel(PointList)
    error('ShiftArr and PointList must be the same size');
end

%This loop generates a master list of clusters and finds their convex hull

MasterL = [];

for i=1:n
    tmp = ShiftArr{i};
    j=1;
    numtmp = numel(tmp);
    while j<=numtmp
        
        pnts = [PointList{i}.pnts(tmp(j).ptIdData,:),j*ones(tmp(j).numPoints)];
        if tmp(j).numPoints > 2
            [hull,V] = convhull(pnts(:,1:2));
            %stores density per pixel^2
            tmp(j).density=tmp(j).numPoints/V;
            tmp(j).R = sqrt(V/pi); %assumes a circular cluster
            tmp(j).hull = hull;
            j=j+1;
        else
            %if a cluster contains 2 or fewer points it is removed from
            %consideration
            tmp(j)=[];
            numtmp=numtmp-1;
        end
    end
    
MasterL = vertcat(MasterL, tmp);    
end


% Make a distance matrix comparing all the cluster centers
cntrs = vertcat(MasterL.ptClusterCenter);
numcnt = numel(cntrs);
R = vertcat(MasterL.R)/2;
%A rough estimate of overlapping clusters. If the two centers are less than
% the sum of the two 'circular' radi over 2.
R2 = repmat(R,[1,numct])+repmat(R',[numct,1]); 
DM = pdist(cntrs);
DM = DM - R2;
DM(DM >0) = Nan;


i=1;
while i <= cntrs
    %loops over the clusters and determines which to merge by DM < 0
    %Once 
end


end



