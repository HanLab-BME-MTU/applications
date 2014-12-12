%   Script to plot cluster densities and fractions over iterations. The
%   desired cluster statistics struct (clusterStat) needs to be loaded into
%   the workspace.
%
%   NOTE: this function was previously called procScript.
%
%   Robel Yirdaw, 02/12/14
%

%Determine number of iterations
numIters = length(clusterStat.clusterCount(1,:,1));
%Set up figure and axes for two plots, densities and fractions.
figure();
ax1 = gca();
set(ax1,'NextPlot','add');
figure();
ax2 = gca();
set(ax2,'NextPlot','add');
%Reset random number generator so that repeated plots will have the same
%sequence of colors. Each cluster size will have a unique color.
rng(0);

%Plot densities and fractions by cluster size
for sizeIndx=1:clusterStat.largestClustSize
    tempH1 = plot(ax1,1:numIters,clusterStat.clusterFracMean(sizeIndx,:));    
    tempH2 = plot(ax2,1:numIters,clusterStat.clusterDensityMean(sizeIndx,:));    
    
    randColor = rand(3,1);
    set(tempH1,'DisplayName',num2str(sizeIndx));
    set(tempH1,'Color',[randColor(1),randColor(2),randColor(3)]);
    
    set(tempH2,'DisplayName',num2str(sizeIndx));
    set(tempH2,'Color',[randColor(1),randColor(2),randColor(3)]);
    
    clear tempH1 tempH2
end

%Adjust axes font, scales and labels
set(ax1,'FontSize',12);
set(ax2,'FontSize',12);
xlabel(ax1,'Iterations','FontSize',12);
xlabel(ax2,'Iterations','FontSize',12);
ylabel(ax1,'Fraction of Clusters','FontSize',12);
ylabel(ax2,'Cluster Density','FontSize',12);
set(ax1,'YScale','log');
set(ax2,'YScale','log');

clear ax1 ax2 numIters randColor sizeIndx


    