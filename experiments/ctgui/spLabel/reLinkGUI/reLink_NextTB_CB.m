function reLink_NextTB_CB(myH,eventdata,handles)
% callback for next_time buttons in reLinkGUI

%check if button is allowed to be activated
noReLink=get(handles.reLink_resetfusionRB,'Value');
if noReLink
    set(myH,'Value',0);
    return %end evaluation here
end

%check if button is being unchecked
myValue=get(myH,'Value');
if myValue==0 %button is being unchecked
    return %end evaluation here
end

%load valuemaps and handleList
guiH=handles.reLinkGUI;
prev_valuemap=GetUserData(guiH,'prev_valuemap');
curr_valuemap=GetUserData(guiH,'curr_valuemap');
next_valuemap=GetUserData(guiH,'next_valuemap');
tag_string=GetUserData(guiH,'tag_string');
buttonH=GetUserData(guiH,'buttonH');

%only myButton should be active
nextH=buttonH{4};
set(nextH,'Value',0);
set(myH,'Value',1);

%find myButtonNumber
myName=get(myH,'Tag');
myNum=str2num(myName(end)); %works for max. 9 buttons, else you have to find the underscore first

%get data on current time and prev time, then decide
%current time
currH=buttonH{3};
curr_status_c=get(currH,'Value');
curr_status=[curr_status_c{1:end}]';
curr_active=find(curr_status);
prevH=buttonH{1};

if ~isempty(curr_active)
    %change string of my button
    statStr=get(currH(curr_active),'String');
    currStr=statStr(1);
    set(myH,'String',currStr,'TooltipString',currStr);
    %change valuemap
    currVal=curr_valuemap(curr_active,2);
    next_other=find(next_valuemap==currVal);
    next_valuemap(myNum)=currVal;
    if ~isempty(next_other)
        next_valuemap(next_other)=0;
        %change string of button previously linked to curr_active
        set(nextH(next_other),'String','?','TooltipString','?');
    end
elseif ~isempty(prevH)
    prev_status_c=get(prevH,'Value');
    prev_status=[prev_status_c{1:end}]';
    prev_active=find(prev_status);
    if ~isempty(prev_active)
        %change string of my button
        prevStr=get(prevH(prev_active),'String');
        set(myH,'String',prevStr,'TooltipString',prevStr);
        %change valuemap
        preVal=prev_valuemap(prev_active);
        next_other=find(next_valuemap(:)==prevVal);
        next_valuemap(myNum)=preVal;
        if ~isempty(next_other)
            next_valuemap(next_other)=0;
            %change string of button previously linked to curr_active
            set(nextH(next_other),'String','?','TooltipString','?');
        end
    end
end

SetUserData(guiH,next_valuemap,1);
