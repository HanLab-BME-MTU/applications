function iClosest = findClosestWindow(windows,point)
%Shittly little function that returns the index of the closest window on the cell-edge (1st band) to the
%input point. (it's assumed the input point is outside the cell or near the
%cell boundary)
%Hunter, 8/2012

nWin = numel(windows);

allD = nan(nWin,1);
for j = 1:nWin    
    if ~isempty(windows{j}) && ~isempty(windows{j}{1})
        currWin = [windows{j}{1}{:}];   
        allD(j) = min(sqrt((point(1) - currWin(1,:)).^2 + (point(2) - currWin(2,:)).^2));        
    end
end

[~,iClosest] = min(allD);

