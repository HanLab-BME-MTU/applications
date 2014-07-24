function reLink_CurrTB_tag_CB(myH,eventdata,handles)
%% callback for curr_time buttons in reLinkGUI

%check mode
tagMode=get(handles.reLink_connectTagsRB,'Value');

%check if button is being unchecked
myValue=get(myH,'Value');
if myValue==0 %button is being unchecked
    return %end evaluation here
end

%load valuemaps and handleList
guiH=handles.reLinkGUI;
curr_valuemap=GetUserData(guiH,'curr_valuemap');
buttonH=GetUserData(guiH,'buttonH');
myName=get(myH,'Tag');
myStr=get(myH,'String');
myNum=str2num(myName(end)); %works for max. 9 buttons, else you have to find the underscore first


switch tagMode
case 0 %fusion mode
    
    %get active spot button
    spotH=buttonH{2};
    spot_status_c=get(spotH,'Value');
    spot_status=[spot_status_c{1:end}]';
    spot_active=find(spot_status);
    
    if isempty(spot_active)
        set(myH,'Value',0); %do nothing
        return %end evaluation here
    else
        %get spot status
        spotName=get(spotH(spot_active),'String');
        spotNum=str2num(spotName);
        %change myString
        myStr(end)=spotName;
        set(myH,'String',myStr,'TooltipString',myStr);
        %change valuemap
        curr_valuemap(myNum,1)=spotNum;
    end
    
    %save data
    SetUserData(guiH,curr_valuemap,1);
    
case 1 %tag mode

    %load additional data
    prev_valuemap=GetUserData(guiH,'prev_valuemap');
    next_valuemap=GetUserData(guiH,'next_valuemap');
    tag_string=GetUserData(guiH,'tag_string');
    currH=buttonH{3};
    prev_active=[];
    
    %only myButton should be active
    set(currH,'Value',0);
    set(myH,'Value',1);
    
    if ~isempty(prev_valuemap)
        prevH=buttonH{1};
        prev_status_c=get(prevH,'Value');
        prev_status=[prev_status_c{1:end}]';
        prev_active=find(prev_status);
        if ~isempty(prev_active)
            %change string of my button
            prevStr=get(prevH(prev_active),'String');
            preVal=prev_valuemap(prev_active);
            %change string of active button
            myStr(1)=prevStr;
            set(myH,'String',myStr,'TooltipString',myStr);
            %change valuemap
            curr_other=find(curr_valuemap(:,2)==preVal);
            curr_valuemap(myNum,2)=preVal;
            if ~isempty(curr_other)
                curr_valuemap(curr_other,2)=0;
                %change string of button previously assigned to myH
                statStr=get(currH(curr_other),'String');
                statStr(1)='?';
                set(currH(curr_other),'String',statStr,'TooltipString',statStr);
            end
        end
    end
    if ~isempty(next_valuemap)&isempty(prev_active) %only relevant of no prev_active
        nextH=buttonH{4};
        next_status_c=get(nextH,'Value');
        next_status=[next_status_c{1:end}]';
        next_active=find(next_status);
        if ~isempty(next_active)
            %change string of next button
            myStr1=myStr(1);
            set(nextH(next_active),'String',myStr1,'TooltipString',myStr1);
            %change valuemap
            myVal=curr_valuemap(myNum,2);
            next_other=find(next_valuemap==myVal);
            next_valuemap(next_active)=myVal;
            if ~isempty(next_other)
                next_valuemap(next_other)=0;
                %change string of button previously linked to curr_active
                set(nextH(next_other),'String','?','TooltipString','?');
            end
        end
    end
    %save data
    SetUserData(guiH,curr_valuemap,1);
    SetUserData(guiH,next_valuemap,1);
end
