function [Output]=meanShift_Merger_poly(ShiftArr, PointList,Overlap)
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
%           Overlap, the area of Overlap required for megering
%                    cluster
%
% Outputs: 
%           Output, A struct containing the merged clusters
%                   
%                   .com  : center of mass of the new cluster
%
%                   .pnts : coordinates of point in a cluster and the
%                           indenty of which data set it came from
%                           [x,y,Id]
%
%                   .hull : indcies of pnts that form convex hull
%
%                   .area : area of cluster
%
%            .composition : a histogram of # of pnts from each element of
%                           ShiftArr
%
%
% Writen 7/19/2013 by Jeffrey Werbin, Harvard Medical School
%

%a few diagnostic variables
countSingle = 0;

%

n = numel(ShiftArr);

if n ~= numel(PointList)
    error('ShiftArr and PointList must be the same size');
end

%This loop generates a master list of clusters and finds their convex hull

MasterL = [];

ClusAdded=[];

for i=1:n
    tmp = ShiftArr{i};
    
    numtmp = numel(tmp);
    tmpRemove = false(size(tmp));
    for j=1:numtmp
        
        pnts = [PointList{i}.pnts(tmp(j).ptIdData,:),j*ones([tmp(j).numPoints,1])];
        tmp(j).pnts = pnts;
        if tmp(j).numPoints > 2
            [hull,V] = convhull(pnts(:,1:2));
            %stores density per pixel^2
            tmp(j).density=tmp(j).numPoints/V;
            tmp(j).R = sqrt(V/pi); %assumes a circular cluster
            tmp(j).hull = pnts(hull,:);
            ClusAdded = vertcat(ClusAdded,[i,j]);
        else
            %if a cluster contains 2 or fewer points it is removed from
            %consideration
            %Note this may be a problem as it discrads these points
            tmpRemove(j)=true;
            countSingle =countSingle+1;
        end
    end
    tmp(tmpRemove)=[];
MasterL = [MasterL, tmp];    
end


%This keeps track of how many clusters are merged
before =0;
after = 0;

% Make a distance matrix comparing all the cluster centers
cntrs = vertcat(MasterL.ptClusterCenter);
numcnt = numel(cntrs(:,1));
R = vertcat(MasterL.R)/2;
%A rough estimate of overlapping clusters. If the two centers are less than
% the sum of the two 'circular' radi over 2.

%Makes the distance matrix a Global variable so it doesn't get copied in
%every step of the recursive search.

R2 = repmat(R,[1,numcnt])+repmat(R',[numcnt,1]); 
DM = squareform(pdist(cntrs));
DM = DM - R2;

Output =[];


%Keeps track of which clusters have already been visited
ListDone =[];
numcnt
for i = 1:numcnt
    
    %Checks if the ith cluster has already been merged and then skips it
    if find(ListDone == i)
        continue;
    end
    
    %loops over the clusters and determines which to consider merging by DM < 0
     
    temp = DM(i,:);
    
    list = find(temp<0); %finds candidate points

    list(list == i)=[];
    
    list = recursiveClumps(list,[i,list],DM);
    list = [i,list];
    
    %Make the shape vector that feeds in to Polygons_intersection
    % note S must be a [1, n] array
    S=repmat(struct('P',struct('x',[],'y',[],'hole',0)),[1,numel(list)]);
    
    for j =1:numel(list)
        temp = MasterL(list(j));
        S(j).P.x = temp.hull(:,1);
        S(j).P.y = temp.hull(:,2);
    end
   
    [Geo] = Polygons_intersection(S,0,1e-6);
    
    %set a cut off for mergers
    % This number is in pixels.
    cut = Overlap;
    
    ToTest = numel(list)+1:numel(Geo);
    MergedTemp = [];
    
    
    
    for j = ToTest
        flag = false(size(MergedTemp));
        if Geo(j).area >= cut
            indexCurrent = Geo(j).index;
            for k = 1:numel(MergedTemp)
                % If the current overlaping section shares an index with
                % another merged section. They are merged and the old
                % section is flaged for removal
                if ~isempty(intersect(indexCurrent,MergedTemp{k}))
                    flag(k)=true;
                    indexCurrent = union(indexCurrent, MergedTemp{k});
                end
            end
            
            MergedTemp(flag) = [];
            MergedTemp = [MergedTemp,{indexCurrent}];
            
        end
    end
    
    %The following prevent isolated clusters from being left out of the
    %Analysis
    
    Tested =[];
        
    for h=1:numel(MergedTemp)
        Tested = union(Tested,MergedTemp{h});
    end

    Left = setdiff([1:numel(list)],Tested);
    for h=1:numel(Left)
        MergedTemp = [MergedTemp,{h}];
    end

    
    
    before = before + numel(list);
    after = after + numel(MergedTemp);
    
    %Merges the mean shift results of all the clusters in MergedTemp and
    %stores them in Output
    
    for j = 1:numel(MergedTemp)
        tmp = struct('pnts',[],'hull',[],'area',[],'composition',[],'com',[]);
        
        %Gets index of merged clusters
        ind = list(MergedTemp{j});
        
        %Merges the list of points 
        pnts = vertcat(MasterL(ind).pnts);
        
        com = mean(pnts(:,1:2));
        
        %finds convex hull and area
        [hull,area]= convhull(pnts(:,1:2));
        
        % n is the number of elements of the shift array. Here we store the
        % number of points in this merged cluster from each element of the
        % shift array
        
        composition = hist(pnts(:,3),1:n);
        
        tmp.pnts = pnts;
        tmp.hull = hull; %indcies to pnts
        tmp.area = area;
        tmp.composition = composition;
        tmp.com = com;
        
        Output = vertcat(Output,tmp);       
        
    end
    
    
    %Updates the list of clusters already handled
    
    ListDone = [ListDone,list];


end

sum(ListDone)

countSingle
before
after

end


function [list]=recursiveClumps(start,total,DM)
% recursively determines all clusters that are potenially clustered
% start is a list of cluster identifiers and returns all groups that maybe
% clustered
%


    

    list = [start];
    n = numel(start);
    nt = numel(total);
    
    if n == 0
        return
    end

    for i=1:numel(start)
        l = find(DM(start(i),:)<0);
        
        %removes any elements that have already been considered
        l = setdiff(l,total);
        
        total = union(total,l); %updates elements considered
        
        l = recursiveClumps(l,total,DM);
        list = union(list,l);
    end
end