function idlist=recalcIdlist(idlist,time,opt,dataProperties,recoverSlist)
%recalcIdlist recalculates connectivities between selected timeframes (similar to spotID)
%
%SYNOPSIS idlist=recalcIdlist(idlist,time,opt)
%
%INPUT    idlist: idlist as produced by the latest version of spotID
%         time  : vector of length 1 to 2 indicating which spots are to be reconnected 
%                   (careful: colors will change in the rest of idlist, anyway)
%         opt   :optional structure with fields
%           opt.weight   : weight of intensity difference vs. norm of displacement to be used
%                          in the time interval specified by t1,t2 ([0...{.5}...1])
%
%OUTPUT   idlist: updated idlist
%
%c: 03/03 Jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


totalNumOfFrames=size(idlist,2);

%fill in time
if length(time)==2
    timeStart=time(1);
    timeEnd=time(2);
    if timeEnd>totalNumOfFrames
        timeEnd=totalNumOfFrames;
    end
else
    timeStart=time(1);
    timeEnd=totalNumOfFrames;
end

%is there opt?
oldWeight=idlist(1).stats.weight;
if nargin==2|isempty(opt)
    newWeight=oldWeight(timeStart:timeEnd);
else
    newWeight=opt.weight;
end

%is there dataProperties?
if nargin<4|isempty(dataProperties)
    dataProperties=[];
end

%is there recoveredSlist?
if nargin>4 & ~isempty(recoverSlist)
    recover = 1;
else
    recover = 0;
end

%-----------initialize input for spotID---------------

%prepare options structure
%calculate vector of weights
weightVector=ones(totalNumOfFrames,1).*oldWeight;
weightVector(timeStart:timeEnd)=newWeight;
%opt-structure: nly start analysis at timeStart -> still not shorten
%weightVector: slist has the same length as idlist, but starts with lots of
%empties
opt.weight=weightVector;
opt.save=0;
opt.verbose=0;

slist(1:totalNumOfFrames)=struct('amp',[],'xyz',[],'deltAmp',[],'sigma',[],'goodTime',[],'nspots',[]);
%read statistics
slist(1).deltAmp=idlist(1).stats.deltAmp;
slist(1).sigma=idlist(1).stats.sigma;

nspots=zeros(totalNumOfFrames,1);
goodTime=zeros(totalNumOfFrames,1);
CoMList=zeros(totalNumOfFrames,3);

%build spotlist from idlist.linklist
%find first timepoint
t1=timeStart;
if isempty(idlist(t1).linklist)
    %lookfor first good frame
    while isempty(idlist(t1).linklist);
        t1=t1+1;
    end
end

%first timepoint: do spotNumber according to colors
[sLinklist,sortIdx]=sortrows(idlist(t1).linklist,4);
slist(t1).xyz=sLinklist(:,9:11);
slist(t1).amp=sLinklist(:,8);
nspots(t1)=size(idlist(t1).linklist,1);
CoMList(t1,:)=idlist(t1).centroid;
slist(t1).detectQ=idlist(t1).info.detectQ_Pix;
slist(t1).trackQ=idlist(t1).info.trackQ_Pix;
%get chi^2 back from linklist (we cannot trust the old "noise"-field)
if size(idlist(t1).linklist,2)>11
    slist(t1).noise = idlist(t1).linklist(sortIdx,12)';
end

%store spotNumbers of first frame (to have Q consistent)
firstSpN = idlist(t1).linklist(:,2);

%field discontinued: slist(t1).noise=idlist(t1).info.noise;
if isfield(idlist(t1).info,'trackerMessage')
    slist(t1).trackerMessage = idlist(t1).info.trackerMessage;
end

%all other timepoints: reconstruct spot coordinates&amplitudes
for t=t1+1:totalNumOfFrames
    if ~isempty(idlist(t).linklist)
        nspots(t)=max(idlist(t).linklist(:,2));
        
        detQ = [];
        nse = [];
        traQ = [];
        for i=1:nspots(t) %readout each spot separately
            rowIdx=find(idlist(t).linklist(:,2)==i);
            slist(t).xyz(i,:)=idlist(t).linklist(rowIdx(1),9:11);
            slist(t).amp(i)=sum(idlist(t).linklist(rowIdx,8));
            CoMList(t,:)=idlist(t).centroid;
            
            %restore QMatrix (sorted by spot)
            detQ = blkdiag(detQ,idlist(t).info.detectQ_Pix( (rowIdx(1)-1)*3+1:rowIdx(1)*3,(rowIdx(1)-1)*3+1:rowIdx(1)*3 ) );
            if ~isempty(idlist(t).info.trackQ_Pix)
                traQ = blkdiag(traQ,idlist(t).info.trackQ_Pix( (rowIdx(1)-1)*3+1:rowIdx(1)*3,(rowIdx(1)-1)*3+1:rowIdx(1)*3 ) );
            end
            
            %get chi^2 back from linklist (we cannot trust the old "noise"-field)
            if size(idlist(t).linklist,2)>11
                nse = [nse,idlist(t).linklist(rowIdx(1),12)];
            end
        end
        
        slist(t).detectQ=detQ;
        slist(t).noise=nse;        
        slist(t).trackQ=traQ;
        if isfield(idlist(t).info,'trackerMessage')
            slist(t).trackerMessage = idlist(t).info.trackerMessage;
        end
        
    else
        %nspots stays zero
        slist(t).amp=[];
    end
end



if recover
    recoverTime = recoverSlist.time;
    slist(recoverTime).amp = recoverSlist.amp;
    slist(recoverTime).xyz = recoverSlist.xyz;
    slist(recoverTime).detectQ = recoverSlist.detectQ;
    slist(recoverTime).trackQ = recoverSlist.trackQ;
    if size(idlist(t1).linklist,2)>11 %backwardCompatibility
        slist(recoverTime).noise = recoverSlist.noise;
    else
        slist(recoverTime).noise = [];
    end
    nspots(recoverTime) = recoverSlist.nspots;
    CoMList(recoverTime,:) = recoverSlist.CoM;
end

goodTime(find(nspots~=0))=1;
slist(1).goodTime=goodTime;
slist(1).CoMList=CoMList;
slist(1).nspots=nspots;

%--------------run spotID-------------------
idlistNew=spotID(slist,opt,dataProperties);

%--------------merge idlists----------------
%

%idlist(timeStart).linklist is taken from old-linklist; linkup from and linkdown to second linklist-new have to be adjusted
%find 2nd time
t2=t1+1; 
lengthIdlist = length(idlist);
%lookfor first good frame
while isempty(idlistNew(t2).linklist)&(t2<totalNumOfFrames+1)
    t2=t2+1;
    if t2>lengthIdlist
        % there is no good frame left! - do not change idlist
        % this might be faster if it is further above in the program, but
        % it is clearer if it's here.
        disp('did not recalc idlist - there was no good frame left')
        return
    end
end

%use of field spot/flag is discontinued - if exist, remove in old idlist
if isfield(idlist,'spot')
    idlist=rmfield(idlist,'spot');
end
if isfield(idlist,'flag')
    idlist=rmfield(idlist,'flag');
end

%if t1 is the first nonempty frame: take whole new idlist, but make sure there
%is the right number of spots in frame 1
firstT = 1;
while isempty(idlist(firstT).linklist)
    firstT = firstT + 1;
end

switch firstT == t1
    case 0 %merge idlists the old way
        
        %if new tag is added (should not happen): stop execution
        if size(idlistNew(t2).linklist,1)~=size(idlist(t1).linklist,1)
            errordlg('error in merging idlists: try recalc from beginning or remove some frames before recalc''ing');
        end
        
        %set linkup to idlist(timeStart).linklist - linklists are sorted now!
        idlistNew(t2).linklist(:,6)=idlist(t1).linklist(:,2);
        
        %set linkdown from idlist(timeStart).linklist
        idlist(t1).linklist(:,7)=idlistNew(t2).linklist(:,2);
        
        %merge idlists; t1 is the last frame we want to keep from old idlist
        idlist=[idlist(1:t1),idlistNew(t1+1:end)];
        
        
    case 1 %make idlistNew(t1) into idlist(t1)
        
        %neat idea, but does not conserve the order of Q
        %coord=idlistNew(t1).linklist(:,9:11);
        %[dummy,dummy,col2] = unique(coord(:,1).*coord(:,2).*coord(:,3)); %multiply to make sure it's the same coord
        
        %set the spotnumbers in first entry of idlist back to the old ones.
        %Careful! if recover, there could be more tags now!
        if recover
            %time1 is 1 or the first good timepoint. So don't worry about
            %Q-consistency. Worry about tag numbers, though
            if size(idlistNew(t1).linklist,1)==length(firstSpN)
                col2 = firstSpN;
            else
                col2 = idlistNew(t1).linklist(:,2);
            end
        else
            col2 = firstSpN;
        end
        %write spotnumbers
        idlistNew(t1).linklist(:,2) = col2;
        
        %linkup (linkdown is already correct)
        idlistNew(t2).linklist(:,6) = col2;
        
        %overwrite idlist, but keep old stats
        stats = idlist(1).stats;
        idlist = idlistNew;
        idlist(1).stats = stats;
end
    
%adjust centroid of timeStart
idlist(t1).centroid=mean(idlist(t1).linklist(:,9:11),1);

%---------make sure that there are no unnecessary colors
linkAll = cat(1,idlist.linklist);

%get list of colors in descending order
allColors = sLinklist(:,4);
allColors= allColors(end:-1:1);

for color = allColors' %e.g. [4,2,1] (for assigns the value of col-i in iteration i)
    %find colors that never appear as isolated tags
    if ~any(linkAll(:,3) == color)
        delColor = color;
        t1 = 1;
        
        %delete entry in list of labels
        labelIdx = bsum2bvec(delColor);
        idlist(1).stats.labelcolor(labelIdx)=[];
        
        %delete color at all timepoints
        for t = 1:size(idlist,2)
            if ~isempty(idlist(t).linklist)
                %find row to delete
                delRow = find(idlist(t).linklist(:,4)==delColor);
                %spotNumber of color2delete
                delSpotNum = idlist(t).linklist(delRow,2);
                %intensity of color2delete
                delInt = idlist(t).linklist(delRow,8);
                %delete row
                idlist(t).linklist(delRow,:) = [];
                
                %adjust all colors
                biggerColIdx = find(idlist(t).linklist(:,4)>delColor);
                idlist(t).linklist(biggerColIdx,4) = idlist(t).linklist(biggerColIdx,4)/2;
                
                %find if this deletes the whole spot (should not happen at all!)
                spotRowIdx = find(idlist(t).linklist(:,2)==delSpotNum);
                if isempty(spotRowIdx) %only tag in spot
                    error('recalcIdlist tried to delete a color that was necessary')
                else %adjust intensity
                    idlist(t).linklist(spotRowIdx,3) = sum(idlist(t).linklist(spotRowIdx,4))*ones(size(spotRowIdx));
                    %intensity is distributed according to remainingInt/sumOfRemainingInts
                    sumInt = sum(idlist(t).linklist(spotRowIdx,8));
                    idlist(t).linklist(spotRowIdx,8) = idlist(t).linklist(spotRowIdx,8)*(1+delInt/sumInt);
                end
                
                %redo spot color
                for i = 1:size(idlist(t).linklist,1)
                    sameIdx = find(idlist(t).linklist(:,2)==idlist(t).linklist(i,2));
                    idlist(t).linklist(i,3) = sum(idlist(t).linklist(sameIdx,4));
                end
                
                %adjust centroid
                idlist(t).centroid = mean(idlist(t).linklist(:,9:11),1);
                
                %sort idlist
                idlist(t).linklist = sortrows(idlist(t).linklist,4);
                
                t1 = t;
            end %if ~isempty(idlist(t))
        end %delete-loop
    end %if
end %for
%-----------------------------------------------------

%update statistics
idlist(1).stats.weight=weightVector;
idlist(1).stats.status{end+1}=[date,': recalc''ed idlist from ',num2str(timeStart)];
idlist(1).stats.maxColor=2*max(idlist(t1).linklist(:,4)); %equivalent to sum(LL(:,4))+1
%if a new tag appeared: update labelcolor
oldLabelColor = idlist(1).stats.labelcolor;
newLabelColor = idlistNew(1).stats.labelcolor;
newLabelColor(1:length(oldLabelColor)) = oldLabelColor;
idlist(1).stats.labelcolor = newLabelColor;
