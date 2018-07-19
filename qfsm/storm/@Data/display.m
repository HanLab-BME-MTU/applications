function display(obj)
% Display selected object properties
disp('Data: =====================');
if ~isempty(obj.roiPosition)
    fprintf('| roiPosition: %.0f %.0f %.0f\n', obj.roiPosition);
end
if ~isempty(obj.roiSize)
    fprintf('| roiSize: %.0f %.0f %.0f\n', obj.roiSize);
end
fprintf('| nPoints: %u\n', obj.nPoints);
fprintf('| nClusters: %u\n', obj.nClusters);
fprintf('| runTime: %s\n', secs2hms(obj.runTime));
disp('===========================');
end



