function reLink_CurrTB_spot_CB(myH,eventdata,handles)
%Callback for spot_TB to select active spot for fusion change

%check if button is allowed to be activated
noFuse=get(handles.reLink_connectTagsRB,'Value');
if noFuse
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
curr_valuemap=GetUserData(guiH,'curr_valuemap');
buttonH=GetUserData(guiH,'buttonH');

%only myButton should be active
set(buttonH{2},'Value',0);
set(myH,'Value',1);

%only tags belonging to this spot should be active
tagH=buttonH{3};
set(tagH,'Value',0);
myName=get(myH,'Tag');
myNum=str2num(myName(end)); %works for max. 9 buttons, else you have to find the underscore first
myTags=find(curr_valuemap(:,1)==myNum);
if ~isempty(myTags)
    set(tagH(myTags),'Value',1);
end
    