function idlist=autoLabel(idlist, intList)
% autolabel tries to label the tags stored in idlist based on intensities.
%
%SYNOPSIS idlist=autoLabel(idlist, intList(optional))
%
%INPUT idlist: list of tag attributes detected by spotdetect and linked by
%               spotID
%      intList: (optional) list of tag intensities
%
%OUTPUT idlist with tags labeled
%
%c: 1/03 Jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%test input
if nargin==1
    %create intList from idlist
    
    %init
    nColors=size(idlist(1).linklist,1);
    nFrames=size(idlist,2);
    intList=zeros(nFrames,nColors);
    
    for i=1:nColors
        for t=1:nFrames
            linklist=idlist(t).linklist;
            if ~isempty(linklist)
                llidx=find(linklist(:,4)==2^(i-1));
                intList(t,i)=linklist(llidx,8);
            end
        end
    end
else
    [nFrames,nColors]=size(intList);
end



%assign labels
switch (nColors<3)
    case 1
        idlist(1).stats.labelcolor(1:nColors)=idlist(1).stats.labellist(1);
    case 0
        
%for fitting: make intList w/o zeros and frame number in the first col
zeroIdx=find(sum(intList,2)==0);
fitIntList=[[1:nFrames]',intList];
fitIntList(zeroIdx,:)=[];

%prepare data for fitting
parameters.xdata=fitIntList(:,1);
for i=1:nColors
    parameters.ydata=fitIntList(:,i+1);
    u0=[parameters.xdata,ones(length(parameters.xdata),1)]\parameters.ydata;
    [u(:,i),sigmaLMS(i)]=leastMedianSquare('(ydata-(u(1)*xdata+u(2))).^2',u0,parameters);
end

u=u';
%"mean" intensities
meanList=u(:,1)*0.5*nFrames+u(:,2);

[dummy, rank]=sort(meanList);
%sort sorts in ascending order!
rank=rank(end:-1:1);

%first two: spb
idlist(1).stats.labelcolor(rank(1:2))=idlist(1).stats.labellist(1);
%all others: cen
idlist(1).stats.labelcolor(rank(3:end))=idlist(1).stats.labellist(2);

end

%update history
idlist(1).stats.status{length(idlist(1).stats.status)+1,1}=[date,': autoLabel'];

%add intensity statistics
idlist(1).stats.autoLabel.u=u;
idlist(1).stats.autoLabel.sigmaLMS=sigmaLMS';