function [ medianFilt,outlierIdx,TSFig] = gcaFindOutliersFromMedFilt( dataIn,medFiltWindSize,k,plotFig)
% gcaFindOutliersFromMedFilt
%
% dataIn: the data to be filtered.
% medFiltWindSize: the window size of the median filter
% k: detect outliers.
% OUTPUT:
% medianFilt: the data filtered :
%
% outlierIdx: the index of the outliers
%
% resOutliers: potentially to be fed into the cost function for linking a
% body piece
%k = 3;

if nargin<4 
    plotFig = true; 
end 

% run the typical 1 d median filter and get smoothed time series
medianFilt = medfilt1(dataIn,medFiltWindSize);
% find those points that have significant residuals from this smoothed
% plot
residuals = medianFilt-dataIn;
if plotFig == true
% figure;
TSFig(1).h = setAxis('off');
TSFig(1).name = 'NeuriteOutgrowth_MedFiltOutlierDetect';

scatter(1:length(dataIn),dataIn,'k','filled');
hold on
plot(medianFilt,'Color','k');
end 
% Now it is just determining the appropriate cut-offs for the residuals
% to be considered "outliers"
% can either make this a hard cut-off or use the framework of detectOuliers
% or even just put these residuals as some weight into the cost of
% a body connection upon revisiting if large more probable connection..
%square the residuals
res2 = residuals .^ 2;

%calculate the median of the squared residuals
medRes2 = nanmedian(res2);

%define parameter to remove outliers (see Rousseeuw 1987, p. 202)
magicNumber2 = 1/norminv(.75)^2;

%calculate test-statistic values
testValue = res2 / (magicNumber2 * medRes2);

%determine which observations are inliers and which are outliers
inlierIdx = find(testValue <= k^2);
outlierIdx = find(testValue > k^2);

%Make sure these damn things are row vectors - somehow it returns them
%differently with different input
inlierIdx = inlierIdx(:)';
outlierIdx = outlierIdx(:)';
if plotFig 
    hold on
    scatter(outlierIdx,dataIn(outlierIdx),'r','filled')
    xlabel('Frame Number','FontSize',12,'FontName','Arial');
    ylabel('Neurite Length Measurement (um)', 'FontSize',12,'FontName','Arial');
    set(gca,'FontSize',10,'FontName','Arial');
    
    warning('off','MATLAB:legend:IgnoringExtraEntries')
    legend('Raw Measurements',['Median Filter: (Window Size ' num2str(medFiltWindSize) ')'],...
        'Outlier Detection','box','off', 'Location','BestOutside');
    warning('on','MATLAB:legend:IgnoringExtraEntries')
    
    title('Veil/Stem Outlier Detection' );
    axis([0 length(dataIn) -25 25]);
else
    TSFig = [];
end

end

