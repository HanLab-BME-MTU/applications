function [weightedMean,weightedStd]=brvelocityAnalyze(vel,tps,sigma,step,name,dispVal)
%BRVELOCITYANALYZE return the mean values and plot the distribution 
%Input
%      vel    :  vect, velocity dsitribution
%      tps    :  vect ,time that the velocity reamin in smae state
%      sigma  :  vect ,standart deviation of the velocity at a given time point
%      step   :  int, size of the discrezation for the distribution
%      name   :  str, name for the graph
%Output
%      weightedMean : mean value of vel. wieghted by the time
%      weightedStd  : std of the mean value of vel. wieghted by the time
discret=[0:step:max(vel)+step];

nvel = [];
indexG=find(vel>0);
indexS=find(vel<0);

sigMoy=zeros(length(discret)-1,1);
for i=1:length(discret)-1
    if length((find(vel>discret(i)&vel<discret(i+1))))>=2
        nvel(i)=length(find(vel>discret(i)&vel<discret(i+1)));
        
        [meanStd] = mean(sigma(find(vel>discret(i)&vel<discret(i+1))));
        sigMoy(i)=meanStd;%/length(find(vel>discret(i)&vel<discret(i+1)));
    elseif length((find(vel>discret(i)&vel<discret(i+1))))==1 
        a=find(vel>discret(i)&vel<discret(i+1));
        sigMoy(i)=sigma(a);
        nvel(i)=1;
    elseif isempty(find(vel>discret(i)&vel<discret(i+1)))
        sigMoy(i)=0;
        nvel(i)=0;
    end
end

discretPlot=[step/2:step:max(vel)+step/2];
figure('Name','Statistic analysze');
[AX,H1,H2]=plotyy(discretPlot,nvel,discretPlot,sigMoy);
set(get(AX(1),'Ylabel'),'string','number of frame');
set(get(AX(2),'Ylabel'),'string','Standart deviation');
xlabel('velocity in \mums^{-1}');
title(['Statistic plot with ' num2str(length(vel)) ' data points, of ' name ]);

set(H1,'LineStyle','-');
set(H2,'LineStyle',':');

    

weightedStd  = weightedStats(sigma, tps,'w');
weightedMean = weightedStats(vel, tps,'w');




 