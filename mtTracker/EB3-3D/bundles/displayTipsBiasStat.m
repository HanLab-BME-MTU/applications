function displayTipsBiasStat(kinTracksOrCell,nameOrCell,varargin)
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
    H=setupFigure(1,1,1,'AspectRatio',1,'AxesWidth',4);
end
% bias should be counted at the moment the MT disappear
if(~iscell(kinTracksOrCell))
    kinTracksOrCell={kinTracksOrCell};
end

tipsCountPerTP=cell(1,length(kinTracksOrCell));
tipsCountPerKin=cell(1,length(kinTracksOrCell));
tipsDistBias=cell(1,length(kinTracksOrCell));
tipsDistBiasP1=cell(1,length(kinTracksOrCell));
tipsDistBiasP2=cell(1,length(kinTracksOrCell));
tipsTime=cell(1,length(kinTracksOrCell));
kinCount=cell(1,length(kinTracksOrCell));
kinIndices=cell(1,length(kinTracksOrCell));
timeCell=cell(1,length(kinTracksOrCell));


for ktIdx=1:length(kinTracksOrCell)
    kinTracks=kinTracksOrCell{ktIdx};
    kinCount{ktIdx}=zeros(1,kinTracks.numTimePoints());

    for kinIdx=1:length(kinTracks)
        kinTrack=kinTracks(kinIdx);
        
        TipsTimes=arrayfun(@(mt) kinTrack.associatedTipsP1(mt).f(kinTrack.associatedTipsP1Idx(mt)), ... 
                                    1:length(kinTrack.associatedTipsP1));
        tipsDistBias{ktIdx}=[tipsDistBias{ktIdx} kinTrack.distTipsP1];
        tipsDistBiasP1{ktIdx}=[tipsDistBiasP1{ktIdx} kinTrack.distTipsP1];
        kinIndices{ktIdx}=[kinIndices{ktIdx} kinIdx*ones(size(TipsTimes))];

        tipsTime{ktIdx}=[tipsTime{ktIdx} TipsTimes];
        TipsTimes=arrayfun(@(mt) kinTrack.associatedTipsP2(mt).f(kinTrack.associatedTipsP2Idx(mt)), ... 
                                    1:length(kinTrack.associatedTipsP2));
        tipsDistBias{ktIdx}=[tipsDistBias{ktIdx} kinTrack.distTipsP2];
        tipsDistBiasP2{ktIdx}=[tipsDistBiasP2{ktIdx} kinTrack.distTipsP2];

        tipsTime{ktIdx}=[tipsTime{ktIdx} TipsTimes];      
        kinCount{ktIdx}(kinTrack.f)=kinCount{ktIdx}(kinTrack.f)+1;
        kinIndices{ktIdx}=[kinIndices{ktIdx} kinIdx*ones(size(TipsTimes))];


    end
    [~,~,~,timeCell{ktIdx},tipsCountPerTP{ktIdx}]=statPerIndx(tipsDistBias{ktIdx},tipsTime{ktIdx},max(tipsTime{1}));
    [~,~,~,~,tipsCountPerKin{ktIdx}]=statPerIndx(tipsDistBias{ktIdx},kinIndices{ktIdx},max(kinIndices{1}));

end
%%
% groupID=arrayfun(@(i) i*ones(size(mtDisappTime{i})),1:length(mtDisappTime),'unif',0);
% gscatter([mtDisappTime{:}], [mtDisappBias{:}],[groupID{:}]);
%%
plot(H(1),vertcat(timeCell{:})',vertcat(tipsCountPerTP{:})');
legend(H(1),nameOrCell);
H=setupFigure(2,2,3,'AspectRatio',1,'AxesWidth',4);
boxplot(H(1),[(tipsDistBias{:})],horzcat(cell2mat(arrayfun(@(i) i*ones(size(tipsDistBias{i})),1:length(tipsDistBias),'unif',0))),'Labels',{'Kin','Rand'});
boxplot(H(2),[(tipsDistBiasP1{:})],horzcat(cell2mat(arrayfun(@(i) i*ones(size(tipsDistBiasP1{i})),1:length(tipsDistBiasP1),'unif',0))),'Labels',{'Kin','Rand'});
boxplot(H(3),[(tipsDistBiasP2{:})],horzcat(cell2mat(arrayfun(@(i) i*ones(size(tipsDistBiasP2{i})),1:length(tipsDistBiasP2),'unif',0))),'Labels',{'Kin','Rand'});

figure
boxplot([(tipsCountPerKin{:})],horzcat(cell2mat(arrayfun(@(i) i*ones(size(tipsCountPerKin{i})),1:length(tipsCountPerKin),'unif',0))),'Labels',{'Kin','Rand'});
title('Close +Tips count');
end

function [means,meds,stds,orderedIndex,counts] = statPerIndx(values,index,maxIndex)

means=zeros(1,maxIndex);
meds=zeros(1,maxIndex);
stds=zeros(1,maxIndex);
counts=zeros(1,maxIndex);

orderedIndex=1:maxIndex;
for t=1:maxIndex
    biasAtTime=(values(index==t));
    means(t)=mean(biasAtTime);
    meds(t)=median(biasAtTime);
    stds(t)=std(biasAtTime);
    counts(t)=sum(index==t);
end
end


