function reLink_PrevTB_CB(myH,eventdata,handles)
% callback for prev_time buttons in reLinkGUI

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
set(buttonH{1},'Value',0);
set(myH,'Value',1);

%find myButtonNumber
myName=get(myH,'Tag');
myNum=str2num(myName(end)); %works for max. 9 buttons, else you have to find the underscore first
myVal=prev_valuemap(myNum);

%find other active buttons in current time and next time and set valuemaps

%current time
currH=buttonH{3};
curr_status_c=get(currH,'Value');
curr_status=[curr_status_c{1:end}]';
curr_active=find(curr_status);
if ~isempty(curr_active)
    %change string of active button
    statStr=get(currH(curr_active),'String');
    statStr(1)=tag_string(myVal);
    set(currH(curr_active),'String',statStr,'TooltipString',statStr);
    %change valuemap
    curr_other=find(curr_valuemap(:,2)==myVal);
    curr_valuemap(curr_active,2)=myVal;
    if ~isempty(curr_other)
        curr_valuemap(curr_other,2)=0;
        %change string of button previously assigned to myH
        statStr=get(currH(curr_other),'String');
        statStr(1)='?';
        set(currH(curr_other),'String',statStr,'TooltipString',statStr);
    end
end

%next time
nextH=buttonH{4};
next_status_c=get(nextH,'Value');
if ~isempty(next_status_c)
    next_status=[next_status_c{1:end}]';
    next_active=find(next_status);
    if ~isempty(next_active)
        %change string of active button
        set(nextH(next_active),'String',tag_string(myVal),'TooltipString',tag_string(myVal));
        %change valuemap
        next_other=find(next_valuemap==myVal);
        next_valuemap(next_active)=myVal;
        if ~isempty(next_other)
            next_valuemap(next_other)=0;
            %change string of button previously assigned to myH
            set(nextH(next_other),'String','?','TooltipString','?');
        end
    end
end

SetUserData(guiH,next_valuemap,1);
SetUserData(guiH,curr_valuemap,1);