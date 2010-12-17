function cellArea=calcCellArea(cellMask,hole)
% subtract only holes that are at the interface, holes in the center of
% the cell will be counted!

% seems like as if the innerMask is without holes already! Simply use
% imFill!

% find a cell with a inner hole!


cellArea=sum(cellMask(:));
for holeId=1:length(hole)
        % create a mask from the hole:
        curveHole=hole{holeId}.curve;
        holeMask = false(max(curveHole(:,2)),max(curveHole(:,1)));
        indHole  = sub2ind(size(holeMask),curveHole(:,2), curveHole(:,1));
        holeMask(indHole)=true;
        holeMask= imfill(holeMask,'holes');
        cellMask= imfill(cellMask,'holes');
        
        [rowsHole,colsHole]=size(holeMask);
        [rowsCell,colsCell]=size(cellMask);
        
        maxRows=max(rowsHole,rowsCell);
        maxCols=max(colsHole,colsCell);
        
        % extend both masks such that they are the same size:
        holeMask(maxRows+1,maxCols+1)=false;
        cellMask(maxRows+1,maxCols+1)=0;
        
        % Area of the hole:
        holeArea=sum(holeMask(:));
        
        % Area of the intersection of hole with cell:
        sectMask=cellMask(holeMask);
        sectArea=sum(sectMask(:));
        
        % only if the sectArea is smaller than the whole holeArea, then
        % the hole is at the boundary and the inner part of it must be
        % subtracted:
        
        if sectArea<holeArea
            cellArea=cellArea-sectArea;
        end
end