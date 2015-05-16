function [ output_args ] = gcaAnalysisMakeTestIndividualMovieDistributionPlot(dataMat,varargin)
% gcaAnalysisMakeTestIndividualMovieDistributionPlot 
% small helper function in the GCAAnalysis Group Report to help visualize
% the distributions for a given group of 

%dataMat: an nxm double array wher n (rows) is the number of measurements
%from each movie and m (columns) is the number of movies in the data set 
% (format compatible with the boxplot- padded by NaNs if the 
%  number of observations observed differs among movies)
% names: the names of the movies for the boxplot
% orderForBoxplot: nx1 double array specifying the order of the boxplot- 
%                   as we sometimes want sorted by neurite outgrowth or some other specified parameter. 
%% Input check 
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataMat')
ip.addParameter('sampleNames',[],@(x) iscell(x));
ip.addParameter('axisFont',12,@(x) ischar(x)); 
ip.addParameter('measurementName','Measurement X',@(x) ischar(x)); 
ip.addParameter('orderForBoxplot',[],@(x) isnumeric(x)); 
ip.parse(dataMat,varargin{:});
p = ip.Results;
%%
names = ip.Results.sampleNames;
orderForBoxplot = ip.Results.orderForBoxplot; 
if ~isempty(orderForBoxplot); 
 % sort the data by array 
 dataMat = dataMat(:,orderForBoxPlot) ; 
 names = names(orderForBoxPlot); 
 end 
%% Set Figure
%fsFigure(0.75); 

%% Detect Outlier Distributions 
subplot(2,2,1:2); 
% Put data into cell format
 samplesCell = arrayfun(@(x) dataMat(:,x),1:size(dataMat,2),'uniformoutput',0); % each column is a sample from a specific cell. 
% Take out NaNs
 samplesCell = cellfun(@(x) x(~isnan(x)),samplesCell,'uniformoutput',0); 
 
% Outlier Detect
 outlierIdx = detectEDFOutliers(samplesCell); 
 xlabel(p.measurementName);

%% Use Anderson Darling test to screen for prototypical distributions  
subplot(2,2,3:4); 
% First Make Boxplot of Sample Distributions
     h1=  boxplot(dataMat,'color','k','notch','on','outlierSize',1,'labelorientation','inline','labels',names);
     
     set(h1(:),'Linewidth',1);
     set(gca,'FontName','Arial','FontSize',ip.Results.axisFont); 
     %set(gca,'XTick',[1:numel(names)]); 
     %set(gca,'XTickLabels',names, 'FontSize',14); 
 

 testDists{1} = 'norm'; 
 testDists{2} = 'exp'; 
 testDists{3} = 'logn';
 testDists{4} = 'weibull';
 %w = warning('off', 'MATLAB:adtest:LargeCoefficient');
 %next step is to see what type of distribution...
 testResults = cellfun(@(x) arrayfun(@(i) adtest(dataMat(:,i),'Distribution',x),1:size(dataMat,2)),testDists,'uniformoutput',0); 
 % test Results will be a cell array for each distribution with a 1 if
 % rejects the null and a 0 if it fails to reject the null hypothesis that
 % the two distributions are the same.
 
 % for each cell find the best distribution match if any
 % testResults = zeros indicates there is evidence it is from that given
 % distribution - therefore take inverse for filtering so that 1 = match.
 testMatrix = ~vertcat(testResults{:}); % should be a nx4 matrix where n is the number of cells of ones and zeros
 
 idxOfDistHitsAll = arrayfun(@(i) find(testMatrix(i,:)),1:4,'uniformoutput',0);
 % idxOfDist Hits is a
 plusX = -5:5:10 ;
 %arrayfun(@(i) testMatrix(:,i),
 % put a text box with the distribution over each hit boxplot
 
% for each distribution
for iDist = 1:4
     if ~isempty(idxOfDistHitsAll{iDist})
     idxOfDistHists = idxOfDistHitsAll{iDist};
 
     arrayfun(@(x) text(idxOfDistHists(x),nanmedian(dataMat(:,x)),testDists{iDist},'FontSize',10),1:length(idxOfDistHists));
     end
end

%% HighLight any outlier distributions from 
 %idxN = 1:size(dataMat,2);
 if ~isempty(outlierIdx)   
 %idxLOut = arrayfun(@(x) idxN==outlierIdx(x),1:length(outlierIdx),'uniformoutput',0);
 %idxLOutFinal = sum(vertcat(idxLOut{:}),1)'; 
 %idxIn = setdiff(idxN,outlierIdx);  
 %dataMat(1:size(dataMat,1),idxIn)= NaN;
 % hold on
 % h22 = boxplot(dataMat,'color','r','colorGroup',idxLOutFinal,'notch','on','outlierSize',1,'labelorientation','inline','labels',names);

 set(h1(:,outlierIdx),'color','r');
 set(h1(:,outlierIdx),'Linewidth',5);
 end 
end

