function [weightedMean,weightedStd]=brvelocityAnalyze(vel,sigmaGlobal,sigma,step,name,dispVal)


discret=[0:step:max(vel)+step];

nvel = [];
indexG=find(vel>0);
indexS=find(vel<0);

sigMoy=zeros(length(discret)-1,1);
for i=1:length(discret)-1
    if length((find(vel>discret(i)&vel<discret(i+1))))>=2
        nvel(i)=length(find(vel>discret(i)&vel<discret(i+1)));
        [meanvel,dummy,stdTemp] = weightedStats(vel(find(vel>discret(i)&vel<discret(i+1))), sigmaGlobal(find(vel>discret(i)&vel<discret(i+1))),'w'); 
        [meanStd,dummy,stdTemp] = weightedStats(sigma(find(vel>discret(i)&vel<discret(i+1))), sigmaGlobal(find(vel>discret(i)&vel<discret(i+1))),'w');
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

    

weightedStd  = weightedStats(sigma, sigmaGlobal,'w');
weightedMean = weightedStats(vel, sigmaGlobal,'w');
 