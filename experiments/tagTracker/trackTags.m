function [idlisttrack]=trackTags(mov,idlist,dataProperties)
%TRACKTAGS main function for chromosome tracker
%
% SYNOPSIS  nsl=trackTags(mov, sl)
%
% INPUT mov : raw microscopy data
%            sl:    list of spots from spotfinder
%
% OUTPUT nsl : new spots list after tracker
%   
% c: 8/1/02 dT


% General Init 
sl_start=1;
sl_end=length(idlist);  % number of timepoints
nsl(sl_end)=idlist(end);    %init new output struct
trackerMessage(1:sl_end)=struct('source',[],'iterations',[],'message',[]);

% number of spots list
for t=1:sl_end
    if ~isempty(idlist(t).linklist)
        nsp(t)=max(idlist(t).linklist(:,2)); 
    else
        nsp(t)=0;
    end;
end;




ctidx=0;
msg={'', '',''};

%find source and target frame pairs
trackPairs=trackingStrategy(idlist);

%init waitbar & info box
maxWaitbar=size(trackPairs,1);
waitbarHandle=mywaitbar(0,[],maxWaitbar,'tracking...');
mHandle=myMessageBox(0,'Initializing...','Tracking',[0 50]);

% main loop
while(~isempty(trackPairs))  % no for loop through list, such that the trackPairs could be adjusted on the fly!
    % take tracked coordinates if possible
    tcsl=idlist(trackPairs(1,:));
    if ~isempty(nsl(trackPairs(1,1)).linklist)
        tcsl(1)=nsl(trackPairs(1,1));
    end;
    msg{1}=['Tracking: ' num2str(trackPairs(1,1)) ' -> ' num2str(trackPairs(1,2))];
    [tempnsl, status]=trackFrame(mov(:,:,:,:,trackPairs(1,:)),tcsl,dataProperties);
    
    ctidx=ctidx+1;
    mywaitbar(ctidx/maxWaitbar,waitbarHandle,maxWaitbar);
    msg{2}=['Number of Iterations: ' num2str(status.iterCt)];
    msg{3}=['Status: ' status.msg];
    mHandle=myMessageBox(mHandle,msg);
    %update values
    nsl(trackPairs(1,2))=tempnsl;
    nsp(trackPairs(1,2))=max(tempnsl.linklist(:,2));
    trackerMessage(trackPairs(1,2)).source{end+1} = num2str(trackPairs(1,1));
    trackerMessage(trackPairs(1,2)).iterations{end+1} = num2str(status.iterCt);
    trackerMessage(trackPairs(1,2)).message{end+1} = status.msg;
    trackPairs=trackPairs(2:end,:);
end;


%fill up non-tracked frames if necessary
for t=1:length(nsl)
    if isempty(nsl(t).linklist)
        nsl(t)=idlist(t);
        % and sort
        if ~isempty(nsl(t).linklist)
            nsl(t).linklist=sortrows(nsl(t).linklist,4);
        end;
    end;
    if ~isempty(nsl(t).linklist)
        nsl(t).centroid=mean(nsl(t).linklist(:,9:11),1);
        nsl(t).info.trackerMessage = trackerMessage(t);
    end;
end;

%correct spot values in linklist
t1=min(find(nsp));
for t2=t1+1:sl_end
    if  ~isempty(nsl(t2).linklist)
        %works only if all linklists are sorted according to tag color!
        nsl(t1).linklist(:,7)=nsl(t2).linklist(:,2);
        nsl(t2).linklist(:,6)=nsl(t1).linklist(:,2);
        t1=t2;
    end
end;


%finishup
close(mHandle);
close(waitbarHandle);
idlisttrack=nsl;
save(['idlisttrack-' nowString],'idlisttrack');
