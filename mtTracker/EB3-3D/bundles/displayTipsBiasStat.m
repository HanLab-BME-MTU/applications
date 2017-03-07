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

tipsCount=cell(1,length(kinTracksOrCell));
tipsDistBias=cell(1,length(kinTracksOrCell));
tipsTime=cell(1,length(kinTracksOrCell));
kinCount=cell(1,length(kinTracksOrCell));
timeCell=cell(1,length(kinTracksOrCell));


for ktIdx=1:length(kinTracksOrCell)
    kinTracks=kinTracksOrCell{ktIdx};
    kinCount{ktIdx}=zeros(1,kinTracks.numTimePoints());

    for kinIdx=1:length(kinTracks)
        kinTrack=kinTracks(kinIdx);
        
        TipsTimes=arrayfun(@(mt) kinTrack.associatedTipsP1(mt).f(kinTrack.associatedTipsP1Idx(mt)), ... 
                                    1:length(kinTrack.associatedTipsP1));
        tipsDistBias{ktIdx}=[tipsDistBias{ktIdx} kinTrack.distTipsP1];
        tipsTime{ktIdx}=[tipsTime{ktIdx} TipsTimes];

        TipsTimes=arrayfun(@(mt) kinTrack.associatedTipsP2(mt).f(kinTrack.associatedTipsP2Idx(mt)), ... 
                                    1:length(kinTrack.associatedTipsP2));
        tipsDistBias{ktIdx}=[tipsDistBias{ktIdx} kinTrack.distTipsP2];
        tipsTime{ktIdx}=[tipsTime{ktIdx} TipsTimes];      
    end
    [~,~,~,timeCell{ktIdx},tipsCount{ktIdx}]=statPerIndx(tipsDistBias{ktIdx},tipsTime{ktIdx},max(tipsTime{1}));
end
%%
% groupID=arrayfun(@(i) i*ones(size(mtDisappTime{i})),1:length(mtDisappTime),'unif',0);
% gscatter([mtDisappTime{:}], [mtDisappBias{:}],[groupID{:}]);
%%
plot(H(1),vertcat(timeCell{:})',vertcat(tipsCount{:})');
legend(H(1),nameOrCell);
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


