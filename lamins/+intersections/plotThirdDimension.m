function [ hpt ] = plotThirdDimension( threeD , threeD2)
%plotThirdDimension Adds an impoint that will plot to the current image

if(nargin < 2)
    threeD2 = [];
end

mainFig = gcf;
mainAx = gca;

hpt = impoint;

plotFig = figure;
plotAx = gca;

% if(~isempty(threeD2))
%     plotFig2 = figure;
%     plotAx2 = gca;
% end

addNewPositionCallback(hpt,@doPlot);

figure(mainFig);

    function doPlot(x)
        data1 = joinColumns(threeD(round(x(2)),round(x(1)),:));
        minData1 = min(data1(:));
        data1 = data1 - min(data1(:));
        data1 = data1 ./ max(data1(:));
        figure(plotFig);
        plot(data1);
        if(~isempty(threeD2))
            data2 = joinColumns(threeD2(round(x(2)),round(x(1)),:));
            x = 1:length(data2);
            data2 = data2 - minData1;
            data2 = data2 ./ max(data2(:));
            hold on;
            plot(x(data2 > 0),data2(data2 > 0),'ro');
            hold off;
        end
        figure(mainFig);
    end


end

