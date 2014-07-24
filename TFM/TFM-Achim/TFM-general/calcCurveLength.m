function curveLength=calcCurveLength(curve,dPts,hole)
% Input: 
% curve:    Nx2 Matrix, with the first colum beeing the x-coordinates and
%           the second column beeing the y-coordinate.
% dPts:     steps, to smooth pixel-curves.
% hole:     sometimes there are holes in the cluster that might want to be
%           substracted from the total length.

% In order to calculate the interfacial length, we take only every 10th
% point:
if isempty(dPts)
    dPts=10;
end


curve_sparse=curve(1:dPts:end,:);

% if the last point is skipped, enter it:
if ~compPts(curve(end,:),curve_sparse(end,:))
    curve_sparse(end+1,:)=curve(end,:);
end

curveLength=sum(sqrt(sum((curve_sparse(2:end,:)-curve_sparse(1:end-1,:)).^2,2)));

% If there are holes, then reduce the curveLength by the length of all
% intersections with holes.
if nargin > 2 && ~isempty(hole)
    holeLength=zeros(1,length(hole));
    for holeId=1:length(hole)
        % Intersect the interface with the hole:
        min_x=min(min(hole{holeId}.curve(:,1)),min(curve(:,1)));
        min_y=min(min(hole{holeId}.curve(:,2)),min(curve(:,2)));
        
        curveHole(:,1)=hole{holeId}.curve(:,1)-min_x+1;
        curveHole(:,2)=hole{holeId}.curve(:,2)-min_y+1;
        
        curveIn(:,1)=curve(:,1)-min_x+1;
        curveIn(:,2)=curve(:,2)-min_y+1;        
        
        max_x=max(max(curveHole(:,1)),max(curveIn(:,1)));
        max_y=max(max(curveHole(:,2)),max(curveIn(:,2)));
        BWmask=zeros(max_y,max_x);
        indHole =sub2ind(size(BWmask),curveHole(:,2), curveHole(:,1));
        indCurve=sub2ind(size(BWmask),curveIn(:,2), curveIn(:,1));       
        BWmask(indHole)=1;
        BWmaskHole = imfill(BWmask,'holes');
        checkVec=BWmaskHole(indCurve);
        interSec=curve(checkVec==1,:);
        if ~isempty(interSec)
            holeLength(holeId)=calcCurveLength(interSec,dPts);
        end
        
        clear curveHole curveIn;
    end
    curveLength=curveLength-sum(holeLength);
end



