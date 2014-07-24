function adguiLeg_redrawLegend
%updates legend for adguiLegend

%get handles&data
adguiH = findall(0,'Tag','adgui');
adguiHandles = guidata(adguiH);
lineStyles = adguiHandles.lineStyles; %{'-.','x';'-','d';':','+';'--','*'};
plotFigureH = adguiHandles.plotFigureH;

legendH = findall(0,'Tag','adguiLegend');
legendHandles = guidata(legendH);
selectedFigure = get(legendHandles.adguiLeg_figure_PD,'Value');

%if you really try hard, you can have a valid handle appear twice - hence the (1)
figureH = plotFigureH(selectedFigure(1));
figureHandles = guidata(figureH);

plotData = figureHandles.plotData;
dataH = figureHandles.plotInfo.dataH;

%dataH has entries in every col for every idlist that was used and in every row
%for the different kinds of data that were plotted. For every nonzero entry on
%dataH there is data stored in plotData(rowIdx,colIdx)
[dataHnumProp,dataHData] = find(dataH);

position = get(legendH,'Position');
numEntries = length(dataHData); % number of entries in legend

newHeight = 60*(1+numEntries);
oldHeight = position(4);
%deltah is positive when the figure size has to get smaller
deltaHeight = oldHeight - newHeight;

oldY = position(2);
newY = oldY + deltaHeight;

%change size of legend
set(legendH,'Position',[position(1),newY,position(3),newHeight]);

%change position of pulldown and associated text
pd_positionC = get(legendHandles.pdH,'Position');
pd_position = cell2mat(pd_positionC);
pd_position(:,2) = pd_position(:,2) - [deltaHeight;deltaHeight];
pd_positionC = mat2cell(pd_position,[1,1],4);
set(legendHandles.pdH,{'Position'},pd_positionC);

%change place of axes, text. First: decide whether to delete or add entries
if deltaHeight > 0
    %figure got smaller, some entries have to be deleted
    delete(legendHandles.axH(numEntries+1:end));
    delete(legendHandles.txtH(numEntries+1:end));
    legendHandles.axH = legendHandles.axH(1:numEntries);
    legendHandles.txtH = legendHandles.txtH(1:numEntries);
    
    %shift remaining axes/texts
    at_positionC = get([legendHandles.axH;legendHandles.txtH],'Position');
    at_position = cell2mat(at_positionC);
    at_position(:,2) = at_position(:,2) - deltaHeight*ones(size(at_position,1),1);
    at_positionC = mat2cell(at_position,ones(1,size(at_positionC,1)),4);
    set([legendHandles.axH;legendHandles.txtH],{'Position'},at_positionC);
    
elseif deltaHeight < 0
    %shift old axes/texts
    at_positionC = get([legendHandles.axH;legendHandles.txtH],'Position');
    at_position = cell2mat(at_positionC);
    at_position(:,2) = at_position(:,2) - deltaHeight*ones(size(at_position,1),1);
    at_positionC = mat2cell(at_position,ones(1,size(at_positionC,1)),4);
    set([legendHandles.axH;legendHandles.txtH],{'Position'},at_positionC);
    
    %add uicontrols
    for i = length(legendHandles.axH)+1:numEntries
        %copy properties: axes
        axData = get(legendHandles.axH(i-1));
        axDataPN = fieldnames(axData);
        axDataPV = struct2cell(axData);
        
        %set new tag
        axDataPV{49} = ['adguiLeg_ax',num2str(i)];
        
        %set new position
        axPos = axDataPV{45};
        axPos(2) = axPos(2) - 60; %new entry is added below
        axDataPV{45} = axPos;
        
        %delete fields 'being deleted', 'type', xyzlabel, title, currentpoint,
        %children
        axDataPN([4,16,23,53,54,64,78,91]) = [];
        axDataPV([4,16,23,53,54,64,78,91]) = [];
        
        %set new axes
        legendHandles.axH(i) = axes;
        set(legendHandles.axH(i),'Units','pixels'); %to prevent ugly conversion (position is updated before units)
        set(legendHandles.axH(i),axDataPN,axDataPV');
        
        %copy properties: txt
        txtData = get(legendHandles.txtH(i-1));
        
        txtDataPN = fieldnames(txtData);
        txtDataPV = struct2cell(txtData);
        
        %set new tag
        txtDataPV{33} = ['adguiLeg_txt',num2str(i)];
        
        %set new position
        txtPos = txtDataPV{27};
        txtPos(2) = txtPos(2) - 60; %new entry is added below
        txtDataPV{27} = txtPos;
        
        %delete fields 'being deleted', 'type', 'extent'
        txtDataPN([2,12,35]) = [];
        txtDataPV([2,12,35]) = [];
        
        %set new uicontrol
        legendHandles.txtH(i) = uicontrol('Style','text');
        set(legendHandles.txtH(i),txtDataPN,txtDataPV')
        
    end %for-loop
end %if

%set new legend text/lines (do in a loop, because I don't want to plot data by
%setting the axes' children) 
for i = 1:numEntries
    %set text
    set(legendHandles.txtH(i),'String',plotData(dataHnumProp(i),dataHData(i)).legendText);
    
    %plot data. LineStyles = {'-.','x';'-','d';':','+';'--','*'};
    styleNum = rem(dataHnumProp(i),4)+1;
    axes(legendHandles.axH(i));
    cla; %get rid of old data
    
    %draw line
    line([0,1],[.5,.5],'Color',extendedColors(adguiHandles.colorList(mod(dataHData(i)-1,size(adguiHandles.colorList,1))+1,:)),...
        'LineStyle',lineStyles{styleNum,1});
    plot(0.5,0.5,'Marker',lineStyles{styleNum,2},'Color',extendedColors(adguiHandles.colorList(mod(dataHData(i)-1,size(adguiHandles.colorList,1))+1,:)));
    
end

%store data
guidata(legendH,legendHandles);