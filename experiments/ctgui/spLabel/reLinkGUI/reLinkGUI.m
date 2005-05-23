function varargout = reLinkGUI(varargin)
% RELINKGUI Application M-file for reLinkGUI.fig
%    FIG = RELINKGUI launch reLinkGUI GUI.
%    RELINKGUI('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 06-Feb-2004 20:17:14

if nargin == 0  % LAUNCH GUI
    
    % if the gui is already open, the user should not be able to reopen it
    % therefore look for open figure and call update if there was an open
    % figure
    rlH = findall(0,'Tag','reLinkGUI');
    if ~isempty(rlH)
        rlHandles = guidata(rlH);
        reLink_updatePB_Callback(rlHandles.reLink_updatePB, [], rlHandles);
        return
    end
    
    
    fig = openfig(mfilename,'reuse');
    
    % Use system color scheme for figure:
    set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
    
    % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    guidata(fig, handles);
    
    if nargout > 0
        varargout{1} = fig;
    end
    
    %initialize GUI
    initGUI(handles);
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    
    try
        [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
    catch
        disp(lasterr);
    end
    
end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.



% --------------------------------------------------------------------
function varargout = reLink_okPB_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.reLink_okPB.

%check for recalc
myTag = get(h,'Tag');
if strcmp(myTag,'reLink_okPB')
    recalc = 0;
else
    recalc = 1;
end

%load data
guiH=handles.reLinkGUI;
curr_valuemap=GetUserData(guiH,'curr_valuemap');
prev_valuemap=GetUserData(guiH,'prev_valuemap');
next_valuemap=GetUserData(guiH,'next_valuemap');
buttonH=GetUserData(guiH,'buttonH');
allTime=GetUserData(guiH,'allTime');
imgFigureH=GetUserData(openfig('labelgui','reuse'),'currentWindow');
idlist=GetUserData(imgFigureH,'idlist');
dataProperties=GetUserData(imgFigureH,'dataProperties');

%find time
prev_time=allTime(1);
curr_time=allTime(2);
next_time=allTime(3);

%number of tags
nTags=size(idlist(curr_time).linklist,1);

%test if all tags are assigned
test_valuemap=[curr_valuemap,next_valuemap];
if any(test_valuemap(:)==0)
    h=warndlg('At least one tag not assigned (->''?'')!','Warning');
    uiwait(h);
    return %end evaluation here
end

%test if all spots are used
delSpot=[];
i=1;
for i=1:length(buttonH{2})
    if ~any(curr_valuemap(:,1)==i)
        delSpot=[delSpot;i];
    end
end

really = [];
if any(delSpot)
    if length(delSpot)>1
        h=warndlg('Only one spot may be deleted at a time!','Warning');
        uiwait(h);
        return %end evaluation here
    end
    %confirm delete spot
    really=questdlg(['Spot',num2str(delSpot),'is unassigned. Do you really want to continue and delete this spot?'],...
        'Are you sure?','delete spot','no','no');
else %do normal ask
    %confirm user input
    %really=questdlg('Do you really want to change these links?','Are you sure?','yes','yes&recalc','no','yes');
end



switch strcmp(really,'no')+2*(strcmp(really,'yes')|strcmp(really,'delete spot'))
    case 1 %no quit%
        return %end evaluation here
    case 2 %quit but no recalc
        %recalc=0;
    case 0 %quit&recalc
        %recalc=1;
        %no recalc if no next_time
        if isnan(next_time)
            h=warndlg('Current time is last frame => no recalc','Warning');
            uiwait(h);
            recalc=0;
        else
            %set start time for recalc: if no fusion changes start next_time, else start curr_time
            if all(curr_valuemap(:,1)==idlist(curr_time).linklist(:,2))
                recalc_start=next_time;
            else
                recalc_start=curr_time;
            end
        end
end

%update idlist(ct).linklist
curr_idl_old=idlist(curr_time).linklist; %remember old current idlist
next_idl_old=zeros(size(curr_idl_old)); %init remember the other two
prev_idl_old=zeros(size(curr_idl_old));
idlist(curr_time).linklist(:,2)=curr_valuemap(:,1);
%update coordinates
for i=1:size(idlist(curr_time).linklist,1)
    rowIdx=find(curr_valuemap(i,1)==curr_idl_old(:,2));
    idlist(curr_time).linklist(i,9:11)=curr_idl_old(rowIdx(1),9:11);
end
idlist(curr_time).linklist(:,4)=2.^(curr_valuemap(:,2)-1);
%multicolor only added at the end

%if deleteSpot: update spot numbers
if any(delSpot)
    rowidx=find(curr_valuemap(:,1)>delSpot);
    if ~isempty(rowidx)
        %if spot #2 has been deleted, #3 becomes #2
        idlist(curr_time).linklist(rowidx,2)=idlist(curr_time).linklist(rowidx,2)-1;
    end
end

%sort idlist(ct).linklist. 

%Readout Q-matrix and noise first
detQSp = [];
traQSp = [];
nseSp = [];
if isempty(delSpot)
    delNum = 0;
else
    delNum = delSpot;
end
%work with curr_idl_old - the other has already been adjusted!
for i = 1:max(curr_idl_old(:,2));
    if ~any(i == delNum) %do not remember entries of deleted spots
        rowIdx = find(curr_idl_old(:,2)==i);
        detQSp = blkdiag(detQSp,idlist(curr_time).info.detectQ_Pix( (rowIdx(1)-1)*3+1:rowIdx(1)*3,(rowIdx(1)-1)*3+1:rowIdx(1)*3 ));
        if size(curr_idl_old,2)>11 %ensure backwardCompatibility
            nseSp = [nseSp;curr_idl_old(rowIdx(1),12)];
        end
        if ~isempty(idlist(curr_time).info.trackQ_Pix) 
            traQSp = blkdiag(traQSp,idlist(curr_time).info.trackQ_Pix( (rowIdx(1)-1)*3+1:rowIdx(1)*3,(rowIdx(1)-1)*3+1:rowIdx(1)*3 ) );
        end
    end
end

%sort idlist
[idlist(curr_time).linklist,curr_sIdx] = sortrows(idlist(curr_time).linklist,4);

%write back the matrices
detQ = [];
traQ = [];
for j=1:nTags
    %the jth tag now was the ith tag before
    spotN = idlist(curr_time).linklist(j,2);
    detQ = blkdiag(detQ,detQSp( (spotN-1)*3+1:spotN*3,(spotN-1)*3+1:spotN*3 ) );
    if size(idlist(curr_time).linklist,2)>11 %backwardCompatibility
        idlist(curr_time).linklist(j,12) = nseSp(spotN);
    end
    if ~isempty(idlist(curr_time).info.trackQ_Pix) 
        traQ = blkdiag(traQ,traQSp( (spotN-1)*3+1:spotN*3,(spotN-1)*3+1:spotN*3 ) );
    end
end

idlist(curr_time).info.detectQ_Pix=detQ;
idlist(curr_time).info.trackQ_Pix=traQ;


%update idlist(pt).linklist
if ~isnan(prev_time)
    prev_idl_old=idlist(prev_time).linklist; %remember old prev linklist
    %     for i=1:nTags
    %         %update linkup, linkdown
    %         idx=find(prev_valuemap(i)==curr_valuemap(:,2)); %find to which curr_tag the prev_tag links
    %         idlist(prev_time).linklist(i,7)=idlist(curr_time).linklist(idx,2); %prev_linkdown
    %         idlist(curr_time).linklist(idx,6)=idlist(prev_time).linklist(i,2); %curr_linkup
    %     end
    %update linkup, linkdown
    idlist(prev_time).linklist(:,7) = idlist(curr_time).linklist(:,2);
    idlist(curr_time).linklist(:,6) = idlist(prev_time).linklist(:,2);
end

%update idlist(nt).linklist
if ~isnan(next_time)
    %     next_idl_old=idlist(next_time).linklist; %remember old next linklist
    %     oldCol=idlist(next_time).linklist(:,4);
    %     idlist(next_time).linklist(:,4)=2.^(next_valuemap(:)-1);
    %     for i=1:nTags
    %         idx=find(next_valuemap(i)==curr_valuemap(:,2)); %find to which curr_tag the prev_tag links
    %         idlist(next_time).linklist(i,6)=idlist(curr_time).linklist(idx,2); %next_linkup
    %         idlist(curr_time).linklist(idx,7)=idlist(next_time).linklist(i,2); %curr_linkdown
    %     end
    %     %update colors (not links) for all subsequent timesteps
    %     newColList=2.^[0:nTags-1]';
    %     for i=1:nTags
    %         newColList(i,2)=idlist(next_time).linklist(find(oldCol==2^(i-1)),4);
    %     end
    %     for ti=next_time+1:length(idlist)
    %         if ~isempty(idlist(ti).linklist)
    %             for i=1:nTags
    %                 idlist(ti).linklist(i,4)=newColList(find(idlist(ti).linklist(i,4)==newColList(:,1)),2);
    %             end
    %         end
    %     end
    
    %all linklists are sorted according to color. Hence we just have to shuffle
    %the rows in the linklists correctly, and then overwrite the colors with the sorted
    %ones
    [sortedColors,shuffleIdx] = sort(2.^(next_valuemap(:)-1));
    for ti = next_time:length(idlist)
        if ~isempty(idlist(ti).linklist)
            %shuffle rows
            idlist(ti).linklist=idlist(ti).linklist(shuffleIdx,:);
            %reassign colors (don't change linkup/linkdown)
            idlist(ti).linklist(:,4)=sortedColors;
            
            
            %update Q-matrices, but not noise: the noise stays with its
            %spot, but Q has to be sorted according to tag number
            detQ = [];
            nse = [];
            traQ = [];
            for j=1:nTags
                i = shuffleIdx(j); %what was in ith position i before goes to the jth now
                detQ = blkdiag(detQ,idlist(ti).info.detectQ_Pix( (i-1)*3+1:i*3,(i-1)*3+1:i*3 ) );
            end
            if ~isempty(idlist(ti).info.trackQ_Pix) 
                for j=1:nTags
                    i = shuffleIdx(j); %what was in ith position i before goes to the jth now
                    traQ = blkdiag(traQ,idlist(ti).info.trackQ_Pix( (i-1)*3+1:i*3,(i-1)*3+1:i*3 ) );
                end
            end
            idlist(ti).info.detectQ_Pix=detQ;
            idlist(ti).info.trackQ_Pix=traQ;
        end
    end
    
    %change linkup/linkdown for next_time
    idlist(curr_time).linklist(:,7) = idlist(next_time).linklist(:,2);
    idlist(next_time).linklist(:,6) = idlist(curr_time).linklist(:,2);
    
end

%add mulitcolor
for ti=curr_time:size(idlist,2)
    if ~isempty(idlist(ti).linklist)
        for i=1:nTags
            sameIdx=find(idlist(ti).linklist(:,2)==i);
            sameCol=sum(idlist(ti).linklist(sameIdx,4));
            idlist(ti).linklist(sameIdx,3)=sameCol*ones(length(sameIdx),1);
        end
    end
end

%------------------update intensities--------------------------------------
%old current spot intensities
for i=1:max(curr_idl_old(:,2))
    sp(curr_time).amp(i)=sum(curr_idl_old(find(curr_idl_old(:,2)==i),8));
end
if any(delSpot)
    sp(curr_time).amp(delSpot)=[];
end


%test if tag intensity still consistent with spot intensity
for i=1:max(idlist(curr_time).linklist(:,2))
    rowIdx=find(idlist(curr_time).linklist(:,2)==i);
    if ~any(sum(idlist(curr_time).linklist(rowIdx,8))~=sp(curr_time).amp)
        %tag intensities not consistent with spot intensity: set tag-int to 0
        idlist(curr_time).linklist(rowIdx,8)=0;
    end
end

if ~isnan(prev_time)
    %go back in time until only separated tags are found or beginning of movie, update intensities
    %of fused spots; if 0-int tags link to multispots, carry 0-int with you
    t1=curr_time;
    t2=prev_time;
    done=0;
    while ~done      
        if ~isempty(idlist(t2).linklist)
            %store spot intensities of t2 & update intensities of multispots
            for i=1:max(idlist(t2).linklist(:,2))
                rowIdx=find(idlist(t2).linklist(:,2)==i);
                sp(t2).amp(i)=sum(idlist(t2).linklist(rowIdx,8));
                if length(rowIdx)>1
                    %get all intensities for updates of multi-t2 from linkdown to t1
                    allInt=idlist(t1).linklist(rowIdx,8);
                    if any(allInt==0)
                        %if any of the t1-intensities found is zero, we can't know the ratios for this
                        %spot
                        idlist(t2).linklist(rowIdx,8)=0;
                    else %update intensities
                        allIntCorr=sp(t2).amp(i)/sum(allInt)*allInt;
                        for j=1:length(rowIdx)
                            idlist(t2).linklist(rowIdx(j),8)=allIntCorr(j);
                        end
                    end
                end
            end
            
            if size(idlist(t2).linklist,1)==max(idlist(t2).linklist(:,2));
                done=1; %stop looking back if only separated tags
            end
            t1=t2;
        end
        t2=t2-1;
        if t2<1
            done=1; %stop looking back if bof
        end
    end %while-loop
    
    %reverse direction and move forward, except if there are still zeros in frame 1
    %(then jump to curr_time/next_time)
    if ~any(idlist(t1).linklist(:,8)==0)
        t2=t1+1; %t1 is the last good t2
        done=0;
        while ~done
            if ~isempty(idlist(t2).linklist)
                for i=1:max(idlist(t2).linklist(:,2))
                    rowIdx=find(idlist(t2).linklist(:,2)==i);
                    if idlist(t2).linklist(rowIdx(1),8)==0 %if one is zero, then all are
                        %get all intensities for updates of multi-t2 from linkup to t1
                        allInt=idlist(t1).linklist(rowIdx,8);
                        if any(allInt==0)
                            error('can''t do intensity update')
                        else %update intensities
                            allIntCorr=sp(t2).amp(i)/sum(allInt)*allInt;
                            for j=1:length(rowIdx)
                                idlist(t2).linklist(rowIdx(j),8)=allIntCorr(j);
                            end
                        end
                    end
                end
                t1=t2;
            end
            t2=t2+1;
            if t2>curr_time %then t1=curr_time and we can move on to the next loop
                done=1; %stop looking forward if eof
            end
        end %while-loop
    end
end

if ~isnan(next_time)
    %go forward in time until only separated tags are found or end of movie, update intensities
    %of fused spots; if 0-int tags link to multispots, carry 0-int with you
    t1=curr_time;
    t2=next_time;
    done=0;
    while ~done      
        if ~isempty(idlist(t2).linklist)
            %store spot intensities of t2 & update intensities of multispots
            for i=1:max(idlist(t2).linklist(:,2))
                rowIdx=find(idlist(t2).linklist(:,2)==i);
                sp(t2).amp(i)=sum(idlist(t2).linklist(rowIdx,8));
                if length(rowIdx)>1
                    %get all intensities for updates of multi-t2 according to sorted tag colors in
                    %t1
                    allInt=idlist(t1).linklist(rowIdx,8);
                    if any(allInt==0)
                        %if any of the t1-intensities found is zero, we can't know the ratios for this
                        %spot
                        idlist(t2).linklist(rowIdx,8)=0;
                    else %update intensities
                        allIntCorr=sp(t2).amp(i)/sum(allInt)*allInt;
                        for j=1:length(rowIdx)
                            idlist(t2).linklist(rowIdx(j),8)=allIntCorr(j);
                        end
                    end
                end
            end
            
            if size(idlist(t2).linklist,1)==max(idlist(t2).linklist(:,2));
                done=1; %stop looking forward if only separated tags
            end
            t1=t2;
        end
        t2=t2+1;
        if t2>size(idlist,2)
            done=1; %stop looking forward if eof
        end
    end %while-loop
    
    %reverse direction and move backward, except if there are still zeros in frame 1
    %(then jump to curr_time)
    if ~any(idlist(t1).linklist(:,8)==0)
        t2=t1-1; %t1 is the last good t2
        done=0;
        while ~done
            if ~isempty(idlist(t2).linklist)
                for i=1:max(idlist(t2).linklist(:,2))
                    rowIdx=find(idlist(t2).linklist(:,2)==i);
                    if idlist(t2).linklist(rowIdx(1),8)==0 %if one is zero, then all are
                        %get all intensities for updates of multi-t2 from linkdown to t1
                        allInt=idlist(t1).linklist(rowIdx,8);
                        if any(allInt==0)
                            error('can''t do intensity update')
                        else %update intensities
                            allIntCorr=sp(t2).amp(i)/sum(allInt)*allInt;
                            for j=1:length(rowIdx)
                                idlist(t2).linklist(rowIdx(j),8)=allIntCorr(j);
                            end
                        end
                    end
                end
                if size(idlist(t2).linklist,1)==max(idlist(t2).linklist(:,2));
                    done=1; %stop looking backwards if only separated tags - we know it's ok from then on
                end
                t1=t2;
            end
            t2=t2-1;
            if t2<1
                done=1; %stop going backward if bof
            end
        end %while-loop
    else
        error('at least one tag has always zero intensity - delete!')
    end
end



%write idlist-status
idlist(1).stats.status{end+1}=[date,': relinked idlist, frame ',num2str(curr_time)];


if recalc
    idlist=recalcIdlist(idlist,recalc_start,[],dataProperties);
else 
    %do nothing, as intensities have already been updated
end

%save idlist
%idlist=projectData(idlist);
SetUserData(imgFigureH,idlist,1);

% delete(handles.reLinkGUI);
labelgui('refresh');


% --------------------------------------------------------------------
function varargout = reLink_cancelPB_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.reLink_cancelPB.


%confirm user input
really=questdlg('Do you really want to quit without saving changes?','Quit?','yes','no','yes');
if strcmp(really,'no')
    return %end evaluation here
end

%remember position
labelGuiH=findall(0,'Tag','labelgui');
positions = GetUserData(labelGuiH,'positions'); %get position struct
positions.relinkPos = get(handles.reLinkGUI,'Position');
SetUserData(labelGuiH,positions,1);

%use delete, not close so that the closreqFcn is not called
delete(handles.reLinkGUI);


% --------------------------------------------------------------------
function varargout = reLink_connectTagsRB_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.reLink_connectTagsRB.

% make sure that if a user clicks a second time, the button is not turned
% off
rbState = get(h,'Value');
if rbState == 0
    set(h,'Value',1);
    return
end

set(handles.reLink_resetfusionRB,'Value',0);

%reset all buttons
buttonH=GetUserData(handles.reLinkGUI,'buttonH');
allH=cat(1,buttonH{1:end});
set(allH,'Value',0);

% --------------------------------------------------------------------
function varargout = reLink_resetfusionRB_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.reLink_resetfusionRB.

% make sure that if a user clicks a second time, the button is not turned
% off
rbState = get(h,'Value');
if rbState == 0
    set(h,'Value',1);
    return
end

set(handles.reLink_connectTagsRB,'Value',0);

%reset all buttons
buttonH=GetUserData(handles.reLinkGUI,'buttonH');
allH=cat(1,buttonH{1:end});
set(allH,'Value',0);


%---------------------------------------------------------------------
function initGUI(reLink_handles)
%initializes GUI: plots buttons, sets UserData etc. 

set(reLink_handles.reLinkGUI,'Visible','off')

%----------get labelgui handles-------------
imgFigureH=GetUserData(openfig('labelgui','reuse'),'currentWindow');
labelGuiH=findall(0,'Tag','labelgui');
timeH=findall(0,'Tag','slider3');

if isempty(imgFigureH)|isempty(labelGuiH)
    delete(reLink_handles.reLinkGUI);
    h=warndlg('No labelgui or no movie loaded','Warning');
    uiwait(h);
    return %end evaluation here
end

%--------------load labelgui data
idlist=GetUserData(imgFigureH,'idlist');
if isempty(idlist)
    delete(reLink_handles.reLinkGUI);
    h=warndlg('No idlist loaded','Warning');
    uiwait(h);
    return %end evaluation here
end
cMap=GetUserData(labelGuiH,'cMap');
cMapFact=size(cMap,1)/idlist(1).stats.maxColor;
curr_time=get(timeH,'Value');

if isempty(idlist(curr_time).linklist)
    delete(reLink_handles.reLinkGUI);
    h=warndlg('Current timepoint contains no data','Warning');
    uiwait(h);
    return %end evaluation here
end   

%--------------get #of spots & #of tags & everything else about curr&prev&next time

%get data on curr_time
curr_linklist=idlist(curr_time).linklist;
curr_nsp=max(idlist(curr_time).linklist(:,2));
curr_ntag=size(curr_linklist,1);
c_colorList=curr_linklist(:,4);
curr_valuemap(:,2)=log2(c_colorList)+1;
curr_valuemap(:,1)=curr_linklist(:,2);
%set curr_time_txt
set(reLink_handles.reLink_currT_txt,'String',['t=',num2str(curr_time)]);
currTH=zeros(curr_ntag,1);
currSH=zeros(curr_nsp,1);

%find #of spots at prev. time
done=0;
step=1;
%prev_nsp=0;
prev_ntag=0;
prev_time=NaN;
prev_valuemap=[];
prev_linklist=[];
while ~done&(curr_time-step>0)
    if ~isempty(idlist(curr_time-step).linklist)
        prev_time=curr_time-step;
        done=1;
        prev_linklist=idlist(prev_time).linklist;
        %prev_nsp=length(idlist(prev_time).spot);
        p_colorList=prev_linklist(:,4);
        prev_ntag=size(prev_linklist,1);
        prev_valuemap=log2(p_colorList)+1;
    end
    step=step+1;
end
set(reLink_handles.reLink_prevT_txt,'String',['t=',num2str(prev_time)]);
prevH=zeros(prev_ntag,1);

%find #of spots at next time
done=0;
step=1;
%next_nsp=0;
next_ntag=0;
next_time=NaN;
next_valuemap=[];
next_linklist=[];
while ~done&(curr_time+step<=length(idlist))
    if ~isempty(idlist(curr_time+step).linklist)
        next_time=curr_time+step;
        done=1;
        next_linklist=idlist(next_time).linklist;
        %next_nsp=length(idlist(next_time).spot);
        n_colorList=next_linklist(:,4);
        next_ntag=size(next_linklist,1);
        next_valuemap=log2(n_colorList)+1;
    end
    step=step+1;
end
set(reLink_handles.reLink_nextT_txt,'String',['t=',num2str(next_time)]);
nextH=zeros(next_ntag,1);

%------------------plot buttons

%distance from left text: 1 char. Distance

%set position parameters for buttons
txtPos = get(reLink_handles.reLink_prevT_txt,'Position');
currFramePos = get(reLink_handles.reLink_currT_frame,'Position');
xend = currFramePos(1)+currFramePos(3)-1;
xnul = txtPos(1)+txtPos(3)+1; %all in chars
xdelta = xend-xnul;
ButtonHeight=txtPos(4);
ButtonWidthTag=(xdelta-curr_ntag+1)/curr_ntag;
ButtonWidthSpot=(xdelta-curr_nsp+1)/curr_nsp;

%initialize strings to print on  buttons
tag_string=char(65:64+curr_ntag); %gives ['A','B','C'...] etc. for colors 1,2,4 etc.

%plot prev_time 
if ~isempty(prev_time)
    prevFramePos = get(reLink_handles.reLink_prevT_frame,'Position');
    ButtonYPos = prevFramePos(2)+(prevFramePos(4)-ButtonHeight)/2;
    for i=1:prev_ntag
        prevButtonColor=cMap(p_colorList(i)*cMapFact,:);
        prevH(i,1)=uicontrol('Style','togglebutton','BackgroundColor',prevButtonColor,...
            'Tag',['PrevTB_',num2str(i)],'Units','characters',...
            'Position',[xnul+(i-1)*(ButtonWidthTag+1),ButtonYPos,ButtonWidthTag,ButtonHeight],...
            'Callback','reLink_PrevTB_CB(gcbo,[],guidata(gcbo))','String',tag_string(prev_valuemap(i)),...
            'TooltipString',tag_string(prev_valuemap(i)),...
            'Parent',reLink_handles.reLinkGUI);
        
        
        
        %if pure color: add strip of mixed color
        multiplicity=length(find(idlist(prev_time).linklist(:,2)==prev_linklist(i,2)));
        if multiplicity>1 
            % read positions in pixel -> change units of button first
            set(prevH(1,1),'Units','pixels')
            buttonPosPix = get(prevH(1,1),'Position');
            % set units back to char
            set(prevH(1,1),'Units','characters');
            % draw strip
            [dummy,dummy,z]=meshgrid(1:buttonPosPix(3),1:round(buttonPosPix(4)/3),1:3);
            stripColor=cMap(round(prev_linklist(i,3)/multiplicity*cMapFact),:);
            strip=stripColor(z);
            set(prevH(i),'CData',strip);
        end
    end
end

ButtonYPos = currFramePos(2)+(currFramePos(4)-ButtonHeight*2)/3*2+ButtonHeight;
%plot curr_time_spot
for i=1:curr_nsp
    
    
    rowIdx=find(curr_linklist(:,2)==i);
    spotColor=curr_linklist(rowIdx(1),3);
    multiplicity=length(rowIdx);
    buttonColor=cMap(round((spotColor)/multiplicity*cMapFact),:);
    currSH(i)=uicontrol('Style','togglebutton','BackgroundColor',buttonColor,...
        'Tag',['CurrTB_spot_',num2str(i)],'Units','characters',...
        'Position',[xnul+(i-1)*(ButtonWidthSpot+1),ButtonYPos,ButtonWidthSpot,ButtonHeight],...
        'Callback','reLink_CurrTB_spot_CB(gcbo,[],guidata(gcbo))',...
        'String',[num2str(i)],...
        'TooltipString',[num2str(i)],...
            'Parent',reLink_handles.reLinkGUI);
    %set units back to char
    set(currSH(i,1),'Units','characters');
end

ButtonYPos = currFramePos(2)+(currFramePos(4)-ButtonHeight*2)/3;
%plot curr_time_tag
for i=1:curr_ntag
    currButtonColor=cMap(c_colorList(i)*cMapFact,:);
    currTH(i)=uicontrol('Style','togglebutton','BackgroundColor',currButtonColor,...
        'Tag',['CurrTB_tag_',num2str(i)],'Units','characters',...
        'Position',[xnul+(i-1)*(ButtonWidthTag+1),ButtonYPos,ButtonWidthTag,ButtonHeight],...
        'Callback','reLink_CurrTB_tag_CB(gcbo,[],guidata(gcbo))',...
        'String',[tag_string(curr_valuemap(i,2)),'-',num2str(curr_valuemap(i,1))],...
        'TooltipString',[tag_string(curr_valuemap(i,2)),'-',num2str(curr_valuemap(i,1))],...
            'Parent',reLink_handles.reLinkGUI);
    
    %if ~pure color: add strip of mixed color
    multiplicity=length(find(idlist(curr_time).linklist(:,2)==curr_linklist(i,2)));
    if multiplicity>1 
        % read positions in pixel -> change units of button first
        set(currTH(1,1),'Units','pixels')
        buttonPosPix = get(currTH(1,1),'Position');
        % set units back to char
        set(currTH(1,1),'Units','characters');
        % draw strip
        [dummy,dummy,z]=meshgrid(1:buttonPosPix(3),1:round(buttonPosPix(4)/3),1:3);
        stripColor=cMap(round(curr_linklist(i,3)/multiplicity*cMapFact),:);
        strip=stripColor(z);
        set(currTH(i),'CData',strip);
    end
end

%plot next_time 
if ~isempty(next_time)
    nextFramePos = get(reLink_handles.reLink_nextT_frame,'Position');
    ButtonYPos = nextFramePos(2)+(prevFramePos(4)-ButtonHeight)/2;
    
    for i=1:next_ntag        
        nextButtonColor=cMap(n_colorList(i)*cMapFact,:);
        nextH(i,1)=uicontrol('Style','togglebutton','BackgroundColor',nextButtonColor,...
            'Tag',['NextTB_',num2str(i)],'Units','characters',...
            'Position',[xnul+(i-1)*(ButtonWidthTag+1),ButtonYPos,ButtonWidthTag,ButtonHeight],...
            'Callback','reLink_NextTB_CB(gcbo,[],guidata(gcbo))','String',tag_string(next_valuemap(i)),...
            'TooltipString',tag_string(next_valuemap(i)),...
            'Parent',reLink_handles.reLinkGUI);
        %set units back to char
        set(nextH(i,1),'Units','characters');
        
        %if pure color: add strip of mixed color
        multiplicity=length(find(idlist(next_time).linklist(:,2)==next_linklist(i,2)));
        %multiplicity=length(idlist(next_time).spot(next_linklist(i,2)).linkup);
        if multiplicity>1 
            % read positions in pixel -> change units of button first
            set(nextH(1,1),'Units','pixels')
            buttonPosPix = get(nextH(1,1),'Position');
            % set units back to char
            set(nextH(1,1),'Units','characters');
            % draw strip
            [dummy,dummy,z]=meshgrid(1:buttonPosPix(3),1:round(buttonPosPix(4)/3),1:3);
            stripColor=cMap(round(next_linklist(i,3)/multiplicity*cMapFact),:);
            strip=stripColor(z);
            set(nextH(i),'CData',strip);
        end
    end
end

%store linkhandles in cell array
buttonH={prevH,currSH,currTH,nextH};

%store time in array
allTime=[prev_time,curr_time,next_time];

%---------if applicable: set position
positions = GetUserData(labelGuiH,'positions'); %get position struct
if isfield(positions,'relinkPos')
    set(reLink_handles.reLinkGUI,'Position',positions.relinkPos)
end


%----------------save user data
SetUserData(reLink_handles.reLinkGUI,prev_valuemap,1);
SetUserData(reLink_handles.reLinkGUI,next_valuemap,1);
SetUserData(reLink_handles.reLinkGUI,curr_valuemap,1);
SetUserData(reLink_handles.reLinkGUI,tag_string,1);
SetUserData(reLink_handles.reLinkGUI,buttonH,1);
SetUserData(reLink_handles.reLinkGUI,allTime,1);


%----------------set view3D
%init time
start_time=prev_time;
if isnan(start_time)
    start_time=curr_time;
end
end_time=next_time;
if isnan(end_time)
    end_time=curr_time;
end

%get view3DH if exist and show data
view3DH=findall(0,'Tag','view3DGUI');
if ~isempty(view3DH)
    view_handles=guidata(view3DH);
    
    %set vie3DH time and mode
    set(view_handles.view3D_tstart_slider,'Value',start_time);
    set(view_handles.view3D_tstart_txt,'String',num2str(start_time));
    set(view_handles.view3D_tend_slider,'Value',end_time)
    set(view_handles.view3D_tend_txt,'String',num2str(end_time));
    set(view_handles.view3D_setTime_PD,'Value',1);
    view3D_refresh;
end


% make gui appear again
set(reLink_handles.reLinkGUI,'Visible','on')
% make reLinkGUI topmost
figure(reLink_handles.reLinkGUI);


% --- Executes on button press in reLink_updatePB.
function reLink_updatePB_Callback(hObject, eventdata, handles)
% hObject    handle to reLink_updatePB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% just restart the gui, but remember the position

% remember position
labelGuiH=findall(0,'Tag','labelgui');
positions = GetUserData(labelGuiH,'positions'); %get position struct
positions.relinkPos = get(handles.reLinkGUI,'Position');
SetUserData(labelGuiH,positions,1);

% restart
delete(handles.reLinkGUI);
reLinkGUI;


