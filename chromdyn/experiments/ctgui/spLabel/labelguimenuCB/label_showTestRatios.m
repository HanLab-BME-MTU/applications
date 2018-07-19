function label_showTestRatios(hObject,eventdata,handles)
%label_showTestRatios is the menu-callback to show the testRatios-figure

%is item checked
isChecked = get(hObject,'Checked');
imgFigureH = GetUserData(handles.labelgui,'currentWindow');

if strcmp(isChecked,'on') %uncheck, close all windows
    set(hObject,'Checked','off');
    testRatioFig = findall(0,'Tag','testRatioFig');
    close(testRatioFig);
    %delete all the saved intFigH
    labelPanelList = findall(0,'Tag','LabelPanel');
    for i = 1:length(labelPanelList)
        SetUserData(labelPanelList(i),[],1,'testRatioFigH');
    end
    return

else %check and open figure
    set(hObject,'Checked','on');
end

%get data
labelguiH = handles.labelgui;

%get movie path
dataFile = GetUserData(imgFigureH,'dataFile');

% get dataProperties
dataProperties = GetUserData(imgFigureH,'dataProperties');

% try to find testRatios-file in dataPath
if ~isempty(dataFile)
    pathName = dataFile.path;
    testRatioFiles = searchFiles('testRatios','',pathName);
    % is there any file with the name testRatios?
    if ~isempty(testRatioFiles)
        % choose good file
        if size(testRatioFiles,1) == 1
            % there is only one file. Easy
            fileName = testRatioFiles{1,1};
        else
            % check whether there is a testRatios-file with the movieName
            fileIdx = strmatch(sprintf('testRatios_%s',dataProperties.name),testRatioFiles(:,1));
            if ~isempty(fileIdx) && length(fileIdx) == 1
                fileName = testRatioFiles{fileIdx,1};
            else
                % let the user choose
                fileIdx = chooseFileGUI(testRatioFiles(:,1));
                fileName = testRatioFiles{fileIdx,1};
            end
        end
    else
        % if no file found - error
        h = errordlg('Sorry, no testRatios-file found');
        uiwait(h);
        % remove all the figures
        set(hObject,'Checked','off');
        testRatioFig = findall(0,'Tag','testRatioFig');
        close(testRatioFig);
        %delete all the saved intFigH
        labelPanelList = findall(0,'Tag','LabelPanel');
        for i = 1:length(labelPanelList)
            SetUserData(labelPanelList(i),[],1,'testRatioFig');
        end
        return
    end
else
    % no dataFile. Still allow to show testRatios
    [fileName,pathName] = uigetfile({'testRatios*','testRatios-file'},'select testRatios file');
end

% load testRatios
load(fullfile(pathName,fileName))

% check that all is ok
if ~exist('testRatios','var') || isempty(testRatios)
    % if no file found - error
        h = errordlg('Sorry, there is a problems with the testRatios-file!');
        uiwait(h);
        % remove all the figures
        set(hObject,'Checked','off');
        testRatioFig = findall(0,'Tag','testRatioFig');
        close(testRatioFig);
        %delete all the saved intFigH
        labelPanelList = findall(0,'Tag','LabelPanel');
        for i = 1:length(labelPanelList)
            SetUserData(labelPanelList(i),[],1,'testRatioFig');
        end
        return
end

% launch figure.
figureHandle = LG_showTestRatios(testRatios, dataProperties);

%set buttonDownFcn, so that a click on the figure changes the current
%timepoint of labelgui
set(gca,'ButtonDownFcn','label_gotoFrame_BDFCN');
%do it also for all lines
lineHCandidates = findall(figureHandle,'Type','line');
set(lineHCandidates,'ButtonDownFcn','label_gotoFrame_BDFCN');

%save handle
SetUserData(imgFigureH,figureHandle,1,'testRatiosFigH');

