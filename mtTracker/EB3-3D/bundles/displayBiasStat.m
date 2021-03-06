function displayBiasStat(kinTracksOrCell,nameOrCell,varargin)
%Plot and compare building 
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('plotHandleArray',[]);
ip.addParameter('placeHolder',[10 35]);
ip.parse(varargin{:});
p=ip.Results;
H=p.plotHandleArray;
if(isempty(H))
    H=setupFigure(2,3,6,'AspectRatio',1,'AxesWidth',4);
end
% bias should be counted at the moment the MT disappear
if(~iscell(kinTracksOrCell))
    kinTracksOrCell={kinTracksOrCell};
end
mtDisappBias=cell(1,length(kinTracksOrCell));
mtDisappDistBias=cell(1,length(kinTracksOrCell));
mtDisappTime=cell(1,length(kinTracksOrCell));
biasAvgCell=cell(1,length(kinTracksOrCell));
biasMedCell=cell(1,length(kinTracksOrCell));
biasStdCell=cell(1,length(kinTracksOrCell));
distAvgCell=cell(1,length(kinTracksOrCell));
distMedCell=cell(1,length(kinTracksOrCell));
distStdCell=cell(1,length(kinTracksOrCell));
distCountCell=cell(1,length(kinTracksOrCell));

kinCount=cell(1,length(kinTracksOrCell));
biasAvgPerKinCell=cell(1,length(kinTracksOrCell));
cutoffDist=400;

timeCell=cell(1,length(kinTracksOrCell));


for ktIdx=1:length(kinTracksOrCell)
    kinTracks=kinTracksOrCell{ktIdx};
    kinCount{ktIdx}=zeros(1,kinTracks.numTimePoints());
    biasAvgPerKinCell{ktIdx}=zeros(1,kinTracks.numTimePoints());

    for kinIdx=1:length(kinTracks)
        kinTrack=kinTracks(kinIdx);
        
        mtDisappEnd=arrayfun(@(mt) mt.f(end), kinTrack.appearingMTP1);
        mtDisappBias{ktIdx}=[mtDisappBias{ktIdx} kinTrack.MTP1Angle'];
        mtDisappDistBias{ktIdx}=[mtDisappDistBias{ktIdx} kinTrack.distP1'];
        mtDisappTime{ktIdx}=[mtDisappTime{ktIdx} mtDisappEnd'];
        if(~isempty(kinTrack.MTP1Angle))
            kinCount{ktIdx}(kinTrack.f)=kinCount{ktIdx}(kinTrack.f)+1;
            biasAvgPerKinCell{ktIdx}(kinTrack.f)=biasAvgPerKinCell{ktIdx}(kinTrack.f)+mean(kinTrack.MTP1Angle);
        end
        
        mtDisappEnd=arrayfun(@(mt) mt.f(end), kinTrack.appearingMTP2);
        mtDisappBias{ktIdx}=[mtDisappBias{ktIdx} kinTrack.MTP2Angle'];
        mtDisappDistBias{ktIdx}=[mtDisappDistBias{ktIdx} kinTrack.distP2'];       
        mtDisappTime{ktIdx}=[mtDisappTime{ktIdx} mtDisappEnd'];
        if(~isempty(kinTrack.MTP2Angle))
            kinCount{ktIdx}(kinTrack.f)=kinCount{ktIdx}(kinTrack.f)+1;
            biasAvgPerKinCell{ktIdx}(kinTrack.f)=biasAvgPerKinCell{ktIdx}(kinTrack.f)+mean(kinTrack.MTP2Angle);
        end
        
    end
    [biasAvgCell{ktIdx},biasMedCell{ktIdx},biasStdCell{ktIdx},timeCell{ktIdx}]=statPerIndx(mtDisappBias{ktIdx},mtDisappTime{ktIdx});
    [distAvgCell{ktIdx},distMedCell{ktIdx},distStdCell{ktIdx},timeCell{ktIdx}]=statPerIndx(mtDisappDistBias{ktIdx},mtDisappTime{ktIdx});
    [~,~,~,~,distCountCell{ktIdx}]=statPerIndx(mtDisappDistBias{ktIdx}(mtDisappDistBias{ktIdx}<cutoffDist),mtDisappTime{ktIdx}(mtDisappDistBias{ktIdx}<cutoffDist));

    biasAvgPerKinCell{ktIdx}=biasAvgPerKinCell{ktIdx}./kinCount{ktIdx};
end
%%
% groupID=arrayfun(@(i) i*ones(size(mtDisappTime{i})),1:length(mtDisappTime),'unif',0);
% gscatter([mtDisappTime{:}], [mtDisappBias{:}],[groupID{:}]);
%%
plot(H(1),vertcat(timeCell{:})',vertcat(biasAvgCell{:})');
title(H(1),'Average angular bias');
legend(H(1),nameOrCell);

plot(H(2),vertcat(timeCell{:})',vertcat(biasMedCell{:})');
title(H(2),'Median angular bias');
legend(H(2),nameOrCell);
%%
% plot(H(3),vertcat(timeCell{:})',vertcat(biasAvgPerKinCell{:})');
% title(H(1),'Kinetochore-averaged ');
% legend(H(3),nameOrCell);
plot(H(4),vertcat(timeCell{:})',vertcat(distAvgCell{:})');
title(H(4),'Average distance bias');

plot(H(5),vertcat(timeCell{:})',vertcat(distMedCell{:})');
title(H(5),'Median distance bias');

plot(H(6),vertcat(timeCell{:})',vertcat(distCountCell{:})');
ylabel(H(6),'Count');
title(H(6),'Disappearance under 400 nm');


end

function [means,meds,stds,orderedIndex,counts] = statPerIndx(values,index)
means=zeros(1,max(index));
meds=zeros(1,max(index));
stds=zeros(1,max(index));
counts=zeros(1,max(index));

orderedIndex=1:max(index);
for t=1:max(index)
    biasAtTime=(values(index==t));
    means(t)=mean(biasAtTime);
    meds(t)=median(biasAtTime);
    stds(t)=std(biasAtTime);
    counts(t)=sum(index==t);
end
end


