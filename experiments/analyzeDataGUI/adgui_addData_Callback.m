function adgui_addData_Callback(hObject, eventdata, handles)
%adds data to existing plot for analyzeDataGUI

%check input
%at least two inputs have to be selected
labelState = get(handles.labelHandles,'Value');
labelState = cat(1,labelState{:});

selectedLabels = find(labelState ~= 1);

if length(selectedLabels) < 2
    h = warndlg('Select at least two axes!','Insufficient Input');
    uiwait(h);
    return
end


necessaryTags = handles.PD_data(labelState,2);
necessaryTags = cat(1,necessaryTags{:});

xyzTagHandles = [handles.xTagHandles(1:necessaryTags(1));...
        handles.yTagHandles(1:necessaryTags(2));...
        handles.zTagHandles(1:necessaryTags(3))];

tagState = [get(xyzTagHandles,'Value')];

if iscell(tagState)
    tagState = cat(1,tagState{:});
end

if any(tagState == 1)
    h = warndlg('Select all necessary tags!','Insufficient Input');
    uiwait(h);
    return
end

%get current figure and its properties
allowedH = cat(1,handles.plotFigureH);
plotFigH = gcf; %will be changed in the future (add PD)

if isempty(plotFigH) | ~ishandle(plotFigH) | ~any(plotFigH == allowedH)
    %lookfor legend and make selected figure active
    legendH = findall(0,'Tag','adguiLegend');
    if ~isempty(legendH)
        legendHandles = guidata(legendH);
        currentFigNum = get(legendHandles.adguiLeg_figure_PD,'Value');
        try
            plotFigH = allowedH(currentFigNum);
        catch
            plotFigH = [];
        end
        if ~ishandle(plotFigH)
            %draw new
            adgui_drawNew_Callback(handles.adgui_drawNew,[],handles);
            return
        end
    else
        %draw new
        adgui_drawNew_Callback(handles.adgui_drawNew,[],handles);
        return
    end
end



%-------------check compatibility

%get data
plotFigHandles = guidata(plotFigH);
plotInfo = plotFigHandles.plotInfo;
plotData = plotFigHandles.plotData;

oldLabelState = plotInfo.labelState;


%get group number from PD_data
compGroupNew = cat(1,handles.PD_data{labelState,3}); 
compGroupOld = cat(1,handles.PD_data{oldLabelState(1,:),3});

%compare
if  ~all(compGroupNew == compGroupOld)
    h = warndlg('You can only add to plot if the data is compatible!','Incompatible Input');
    uiwait(h);
    adgui_drawNew_Callback(handles.adgui_drawNew,[],handles);
    return
end

%-------------------------------




%---------get data to plot
[plotData,currentData,selectedTags] = adgui_calcPlotData(handles,selectedLabels,labelState);
%-------------------------


%plot data

%make figure topmost
figure(plotFigH);
hold on;

%select style: 
% 1) check in plotInfo.dataH which styles are not used yet for the currentData
% 2) check if among these unused styles there are any that have been plotted in
% this figure already
% 3) if the current labelState corresponds to a style found in
% 2, draw with that style; else use first unused style

%get possibleStyles
allowedStyles = [1:size(plotInfo.dataH,1)+1]; %either use an existing style or add a new one

if size(plotInfo.dataH,2) < currentData %nothing has been plotted for this data yet
    %do nothing
else %keep only styles which have not been used yet for currentData
    allowedStyles(find(plotInfo.dataH(:,currentData))) = [];
end


%try to find corresponding style
selectedStyle = allowedStyles(end);
done = 0;
i = 1;
while ~done & i < length(allowedStyles)
    if plotInfo.labelState(allowedStyles(i),:) == labelState'
        selectedStyle = allowedStyles(i);
        done = 1;
    else
        i = i+1;
    end
end

%select style Number
styleNum = rem(selectedStyle,4)+1;

%get line styles
styles = handles.lineStyles; %{'-.','x';'-','d';':','+';'--','*'};

switch length(selectedLabels)
    case 2 %2D-plot
        dataLength = min(length(plotData.xData),length(plotData.yData));
        dataH = plot(plotData.xData(1:dataLength),plotData.yData(1:dataLength),...
            'Color',extendedColors(handles.colorList(mod(currentData-1,size(handles.colorList,1))+1,:)),...
            'LineStyle',styles{styleNum,1},'Marker',styles{styleNum,2});
        xlabel(plotData.x_Label);
        ylabel(plotData.y_Label);
        
        %plot error bars if option selected
        if get(handles.adgui_errorbar_check,'Value') == 1
            xError =  ~isempty(plotData.xStats);
            yError =  ~isempty(plotData.yStats);
            
            switch xError + 2*yError
                %outcomes: 0 -> no bars
                %          1 -> only xBars, write zeros for yBars
                %          2 -> only yBars
                %          3 -> both xBars and yBars
                case 0
                    %no bars
                case 1
                    bars = [plotData.xStats.sigma(1:dataLength), zeros(dataLength,1)];
                    myErrorbar(plotData.xData(1:dataLength),plotData.yData(1:dataLength),bars);
                case 2
                    bars = [plotData.yStats.sigma(1:dataLength)];
                    myErrorbar(plotData.xData(1:dataLength),plotData.yData(1:dataLength),bars);
                case 3
                    bars = [plotData.xStats.sigma(1:dataLength), plotData.yStats.sigma(1:dataLength)];
                    myErrorbar(plotData.xData(1:dataLength),plotData.yData(1:dataLength),bars);
            end
                    
        end
        
    case 3
        dataLength = min([length(plotData.xData),length(plotData.yData),length(plotData.zData)]);
        dataH = plot3(plotData.xData(1:dataLength),plotData.yData(1:dataLength),plotData.zData(1:dataLength),...
            'Color',extendedColors(handles.colorList(currentData,:)),'Marker',styles{styleNum,2});
        xlabel(plotData.x_Label);
        ylabel(plotData.y_Label);
        zlabel(plotData.z_Label);
        set(gca,'XGrid','on','YGrid','on','ZGrid','on');
end

%store info in figure. PlotData is a #_of_plotted_property x #_of_data_in_dataPD array with fields
% labelState (state of PD's)
% selectedLabels (which of the PD's are active)
% selectedTags (3 by 1 cell containing the tag numbers)
% plotData
% legendText (dataName;x/y/z-label)
%
% plotInfo is an array with the fields 
% dataH (dataHandles for the fields in plotData)
% labelState (n by 2|3 vector of labelStates)


plotFigHandles = guidata(plotFigH);

%selectedStyle is the rowIndex for plotInfo and plotData
plotInfo.labelState(selectedStyle,:) = labelState; 
plotInfo.dataH(selectedStyle,currentData) = dataH;
%plotInfo.switchState = switchState;

plotFigHandles.plotData(selectedStyle,currentData).data = plotData;
plotFigHandles.plotData(selectedStyle,currentData).selectedLabels = selectedLabels;
plotFigHandles.plotData(selectedStyle,currentData).legendText = {handles.data(currentData).dataProperties.name;plotData.legendText};
plotFigHandles.plotData(selectedStyle,currentData).selectedTags = selectedTags;

plotFigHandles.plotInfo = plotInfo;

guidata(plotFigH,plotFigHandles);


%if there is a legend being shown: update it
legendH = findall(0,'Tag','adguiLegend');
if ~isempty(legendH)
    legendHandles = guidata(legendH);
    
    currentFig = find(plotFigH == allowedH);
    
    %set string to active figure
    set(legendHandles.adguiLeg_figure_PD,'Value',currentFig);
    
    %update legend
    adguiLeg_redrawLegend;
end