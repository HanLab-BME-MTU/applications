function cellArea=calcCellArea(cellMask,hole)
% subtract only holes that are at the interface, holes in the center of
% the cell will be counted!

% Since the innerMask is contains any holes (meaning that it has zeros at
% hole coordinates), we can simply use imfill to fill up hole in the inner
% of the mask, holes located at the interface will persist and won't be
% filled up. Thus we, account for the area of inner holes but discount
% areas at the interface.

% INPUT: the second input 'hole' is obsolete.

cellArea=sum(sum(imfill(cellMask,'holes')));


% This is a complicated way of what is be done by the one line above
% cellArea=sum(cellMask(:));
% for holeId=1:length(hole)
%         % create a mask from the hole:
%         curveHole=hole{holeId}.curve;
%         holeMask = false(max(curveHole(:,2)),max(curveHole(:,1)));
%         indHole  = sub2ind(size(holeMask),curveHole(:,2), curveHole(:,1));
%         holeMask(indHole)=true;
%         holeMask= imfill(holeMask,'holes');
%         cellMask= imfill(cellMask,'holes');
%         
%         [rowsHole,colsHole]=size(holeMask);
%         [rowsCell,colsCell]=size(cellMask);
%         
%         maxRows=max(rowsHole,rowsCell);
%         maxCols=max(colsHole,colsCell);
%         
%         % extend both masks such that they are the same size:
%         holeMask(maxRows+1,maxCols+1)=false;
%         cellMask(maxRows+1,maxCols+1)=0;
%         
%         % Area of the hole:
%         holeArea=sum(holeMask(:));
%         
%         % Area of the intersection of hole with cell:
%         sectMask=cellMask(holeMask);
%         sectArea=sum(sectMask(:));
%         
%         % only if the sectArea is smaller than the whole holeArea, then
%         % the hole is at the boundary and the inner part of it must be
%         % subtracted. If sectArea==holeArea, then the hole is located in
%         % the inner of the cell and then it has to be added.
%         
%         if sectArea<holeArea
%             cellArea=cellArea-sectArea;
%         elseif sectArea==holeArea
%             cellArea=cellArea+sectArea;
%         end            
% end