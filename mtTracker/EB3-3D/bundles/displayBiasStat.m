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
    H=setupFigure(1,1,1,'AspectRatio',1,'AxesWidth',4);
end
% bias should be counted at the moment the MT disappear
if(~iscell(kinTracksOrCell))
    kinTracksOrCell={kinTracksOrCell};
end
mtDisappBias=cell(1,length(kinTracksOrCell));
mtDisappTime=cell(1,length(kinTracksOrCell));
biasAvgCell=cell(1,length(kinTracksOrCell));
biasStdCell=cell(1,length(kinTracksOrCell));
timeCell=cell(1,length(kinTracksOrCell));


for ktIdx=1:length(kinTracksOrCell)
    kinTracks=kinTracksOrCell{ktIdx};
    for kinIdx=1:length(kinTracks)
        kinTrack=kinTracks(kinIdx);
        
        mtDisappEnd=arrayfun(@(mt) mt.f(end), kinTrack.appearingMTP1);
        mtDisappBias{ktIdx}=[mtDisappBias{ktIdx} kinTrack.MTP1Angle'];
        mtDisappTime{ktIdx}=[mtDisappTime{ktIdx} mtDisappEnd'];
        
        
        mtDisappEnd=arrayfun(@(mt) mt.f(end), kinTrack.appearingMTP2);
        mtDisappBias{ktIdx}=[mtDisappBias{ktIdx} kinTrack.MTP2Angle'];
        mtDisappTime{ktIdx}=[mtDisappTime{ktIdx} mtDisappEnd'];
        
    end
    biasAvg=zeros(1,kinTracks.numTimePoints());
    biasStd=zeros(1,kinTracks.numTimePoints());
    for t=1:kinTracks.numTimePoints()
        biasAtTime=(mtDisappBias{ktIdx}(mtDisappTime{ktIdx}==t));
        biasAvg(t)=mean(biasAtTime);
        biasStd(t)=std(biasAtTime);
    end
    biasAvgCell{ktIdx}=biasAvg;
    biasStdCell{ktIdx}=biasStd;
    timeCell{ktIdx}=1:kinTracks.numTimePoints();
end
%%
groupID=arrayfun(@(i) i*ones(size(mtDisappTime{i})),1:length(mtDisappTime),'unif',0);
gscatter([mtDisappTime{:}], [mtDisappBias{:}],[groupID{:}]);
%%
plot(H(1),vertcat(timeCell{:})',vertcat(biasAvgCell{:})');
legend(H(1),nameOrCell);

end

