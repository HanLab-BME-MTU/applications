function label_gotoFrame_BDFCN
%ButtonDownFunction for intensity plot and displacement plot that lets the user goto a timepoint in labelgui by clicking on the figure%

%-------get all the necessary handles


%get the labelguiHandles
labelguiH = findall(0,'Tag','labelgui');
if isempty(labelguiH)
  return
end
label_handles = guidata(labelguiH); %handle structure of labelGui 

%get the handle of the calling function (should be of type axes or line)
callingH = gcbo;

%to find the right timepoint, we have to get at the dataLine's XData
%hence check for the type of the callingH. Since a line could be both
% data or errorBar, find first the axesHandle
%and then go back to the lineH

callingHType = get(callingH,'Type');

if strcmp(callingHType,'axes')
    axesH = callingH;
else
    %we don't have to worry about other possible types - there are none with that BDFCN
    axesH = get(callingH,'Parent');
end

lineHCandidates = get(axesH,'Children');

%find one good one

goodLines = findall(lineHCandidates,'Type','line','Tag','');
    
if isempty(goodLines)
    error('no valid lines found for current axes!')
end

%take all lines
lineH = goodLines;

%------end get handles




%------find the currentPoint and the closest XData

%extract XData: take all xData
xData = [];
for i = 1:length(lineH)
xData = [xData;get(lineH(i),'XData')'];
end
xData = unique(xData);

%get current point
currentPoint = get(axesH,'CurrentPoint');
%currentX is the first currentPoint - coordinate
currentX = currentPoint(1);

%find minDistance between currentPointX and XData (no problem, since we're using timepoints)
[val, timePointIdx] = min(abs(xData-currentX));

timePoint = xData(timePointIdx);

%------end find timePoint



%------set timePoint in LabelGUI

%set slider value before calling the function
set(label_handles.slider3,'Value',timePoint);
        
%call the timeSlider function - it should do all the refresh calls by itself
labelgui('slider3_Callback',label_handles.slider3,[],label_handles);

%------end set timePoint

