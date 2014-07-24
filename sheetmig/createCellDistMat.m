function [cellDistFromEdge,distGrad]=createCellDistMat(sheetMask,r,dPix,toDoList,sglCell)
% This function calculates the distance matrix for the croped image
% INPUT
%   sheetMask   : The mask of the cell sheet which might be padded with
%                 zeros.
%   r           : The dx to calculate the discrete gradient of the distance
%                 matrix
%   dPix        : The uniform zero padding of sheetMask.

if nargin < 4 || isempty(toDoList);
    toDoList=1:length(sheetMask);
end

if nargin < 5 || isempty(sglCell);
    sglCell=0;
end

cellDistFromEdge(toDoList)= struct('mat', zeros(size(sheetMask(1).mat)));
if nargout>1
    distGrad(toDoList)=struct('xmat', zeros(size(sheetMask(1).mat)),'ymat', zeros(size(sheetMask(1).mat)));
end

% Doesn't help either:
% cellDistFromEdge = repmat(struct('mat', zeros(size(sheetMask(1).mat))), [1 toDoList(end)]);
% if nargout>1
%     distGrad(toDoList)=repmat(struct('xmat', zeros(size(sheetMask(1).mat)),'ymat', zeros(size(sheetMask(1).mat))),[1 toDoList(end)]);
% end

if sglCell
    for frame=toDoList
        cellDistFromEdge(frame).mat=zeros(size(sheetMask(frame).mat));
        distGrad(frame).xmat=zeros(size(sheetMask(frame).mat));
        distGrad(frame).ymat=zeros(size(sheetMask(frame).mat));
    end
else    
    for frame=toDoList
        text=['Processing ',num2str(toDoList(end)),' frames'];
        progressText(frame/toDoList(end),text);
        
        % crop the mask as before:
        [rowsOrg, colsOrg]=size(sheetMask(frame).mat);
        sheetMaskCrop =sheetMask(frame).mat((1+dPix(frame)):(rowsOrg-dPix(frame)),(1+dPix(frame)):(colsOrg-dPix(frame)));
        
        % shift distance Matrix
        D=bwdist(~sheetMaskCrop);
        linIdSmall = (1:numel(D))';
        [rpos,cpos]=ind2sub(size(D),linIdSmall);
        rpos=rpos+dPix(frame);
        cpos=cpos+dPix(frame);
        
        % cellDistFromEdge(frame).mat=zeros(size(sheetMask(frame).mat));
        linInd=sub2ind(size(cellDistFromEdge(frame).mat),rpos,cpos);
        cellDistFromEdge(frame).mat(linInd)=D(linIdSmall);
        
        % Calculate and shift the gradient:
        % Spacing for calculating the gradient:
        
        if nargout>1
            [DX,DY] = gradient(D,r);
            % distGrad(frame).xmat=zeros(size(sheetMask(frame).mat));
            % distGrad(frame).ymat=zeros(size(sheetMask(frame).mat));
            distGrad(frame).xmat(linInd)=DX(linIdSmall);
            distGrad(frame).ymat(linInd)=DY(linIdSmall);
        end
    end
end
