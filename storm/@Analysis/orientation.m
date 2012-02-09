function or = orientation(obj,movingAverageSpan,displayText,maxHistValue,zoom)

% Tangent of all the points
[~,normalT] = cellfun(@tangentBezier, ...
    obj.data.modelBezCP, ...
    obj.data.modelProj, ...
    'UniformOutput',false);

% Orientation of the individual models
or = cellfun(@(a) sum(a(:,1:2),1)/norm(sum(a(:,1:2),1)),normalT,'UniformOutput',false);
or = cell2mat(or);

% Flip vectors if they point to negative y-values
inverseOr = or(:,2) < 0;
or(inverseOr,:) = -or(inverseOr,:);

modelLength = obj.data.modelLength;

% Compute the angle
[theta,~,~] = cart2sph(or(:,1),or(:,2),zeros(size(or,1),1));

% Compute the histogram
nBins = ceil(sqrt(numel(theta)));
binSize = (max(theta)-min(theta))/nBins;
histEdges = min(theta):binSize:max(theta);
histValues = zeros(1,nBins);
for b=1:nBins-1
    binModelIdx = theta>=histEdges(b) & theta<histEdges(b+1);
    histValues(1,b) = sum(modelLength(binModelIdx));
end
b = nBins;
binModelIdx = theta>=histEdges(b) & theta<=histEdges(b+1);
histValues(1,b) = sum(modelLength(binModelIdx));
% binCenters = min(theta)+binSize/2:binSize:max(theta)-binSize/2;
histValues = smooth(histValues,movingAverageSpan)';

% Build plot data
histValues4Plot = reshape([histValues;histValues;zeros(1,numel(histValues))],3*numel(histValues),1);
binCenters4Plot = reshape([histEdges(1:end-1);histEdges(2:end);histEdges(2:end)],3*numel(histValues),1);
histValues4Plot = [0;histValues4Plot]*zoom;
binCenters4Plot = [binCenters4Plot(1);binCenters4Plot];

% Flip the histogram horizontally
histValues4Plot = histValues4Plot(end:-1:1);

% Set lower bound for the maximum value of the histogram
polarCustom([255 255 255]/255,[binCenters4Plot;0],[histValues4Plot;maxHistValue],'b');

% Remove all the lines
delete(findall(gca,'type','line'));

hold on;

% Plot the main histogram
polarCustom([255 255 255]/255,binCenters4Plot,histValues4Plot,'b');

% Plot the point symmetric part of the histogram
polarCustom([255 255 255]/255,binCenters4Plot+pi,histValues4Plot,'b');

% Plot the predominant orientation
[~,idx] = max(histValues(end:-1:1));
h = polarCustom([255 255 255]/255,binCenters4Plot((idx-1)*3+1:(idx-1)*3+4),histValues4Plot((idx-1)*3+1:(idx-1)*3+4),'r');
set(h,'LineWidth',2);

% Plot the point symmetric part of the predomainant orientation
h = polarCustom([255 255 255]/255,binCenters4Plot((idx-1)*3+1:(idx-1)*3+4)+pi,histValues4Plot((idx-1)*3+1:(idx-1)*3+4),'r');
set(h,'LineWidth',2);

% Remove all the text
t = findall(gcf,'type','text');
if ~displayText
    delete(t);
end

hold off;

end