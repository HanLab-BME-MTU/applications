function [stats,h] = plotShortVSLongLivedCPP(lftRes,threshold)
% plotShortVSLongLivedCPP produce a bar plot of short lived vs long lived
% CCP using a user defined threshold in seconds.
% Philippe Roudot, Marcel Mettlen 2016

thresholdBin=find(lftRes.t==threshold);

h=setupFigure(1,1,1,'AspectRatio',1);
counts=[sum(lftRes.lftHistCCP(:,1:thresholdBin),2),sum(lftRes.lftHistCCP(:,thresholdBin:end),2)];
boxplot(h,counts)
names={['CCP <' num2str(threshold) ' s'],['CCP >' num2str(threshold) ' s']};
set(h,'xtick',[1:2],'xticklabel',names)

stats.mean=mean(counts);
stats.median=median(counts);
stats.std=std(counts);
stats.min=min(counts);
stats.max=max(counts);