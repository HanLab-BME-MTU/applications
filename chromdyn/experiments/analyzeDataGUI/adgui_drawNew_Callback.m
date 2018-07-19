function adgui_drawNew_Callback(hObject, eventdata, handles)

%opens a new figure and draws the specified plot in analyzeDataGUI


%------------test "input"

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


%---------get data to plot
[plotData,currentData,selectedTags] = adgui_calcPlotData(handles,selectedLabels,labelState);
%-------------------------


%plot data
plotFigH = figure('Name',handles.data(currentData).dataProperties.name);

switch length(selectedLabels)
    case 2 %2D-plot
        dataLength = min(length(plotData.xData),length(plotData.yData));
        dataH = plot(plotData.xData(1:dataLength),plotData.yData(1:dataLength),...
            'Color',extendedColors(handles.colorList(mod(currentData-1,size(handles.colorList,1))+1,:)),...
            'Marker','d','LineStyle','-');
        xlabel(plotData.x_Label);
        ylabel(plotData.y_Label);
        grid on;
        hold on;
        
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
            'Color',extendedColors(handles.colorList(currentData,:)),'Marker','d');
        xlabel(plotData.x_Label);
        ylabel(plotData.y_Label);
        zlabel(plotData.z_Label);
        grid on;
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

plotInfo.labelState(1,:) = labelState;
plotInfo.dataH(1,currentData) = dataH;
%plotInfo.switchState = switchState;

plotFigHandles.plotData(1,currentData).data = plotData;
plotFigHandles.plotData(1,currentData).selectedLabels = selectedLabels;
plotFigHandles.plotData(1,currentData).legendText = {handles.data(currentData).dataProperties.name;plotData.legendText};
plotFigHandles.plotData(1,currentData).selectedTags = selectedTags;

plotFigHandles.plotInfo = plotInfo;

guidata(plotFigH,plotFigHandles);

%store figureHandle in adgui
numData = size(handles.plotFigureH,2);
handles.plotFigureH(numData+1) = plotFigH;
guidata(hObject,handles);

%if there is a legend beign shown: update it
legendH = findall(0,'Tag','adguiLegend');
if ~isempty(legendH)
    legendHandles = guidata(legendH);
    %get plotH
    allPlotHandles = handles.plotFigureH;
    
    numberTitleStr = {'Figure No. '};
    numFig = length(allPlotHandles);
    
    %create PD-String
    pdString = cellstr([char(numberTitleStr(ones(numFig,1))),num2str([allPlotHandles]')]);
    
    %set string
    set(legendHandles.adguiLeg_figure_PD,'String',pdString,'Value',numData+1);
    
    %update legend
    adguiLeg_redrawLegend;
end

