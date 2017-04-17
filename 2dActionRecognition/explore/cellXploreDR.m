function [cxFig handles] = cellXploreDR(cellDataSet, varargin)
%CELLXPLORER Interactive display of movies associated with point on a 2D scatter plot.
%
% Inputs:
% 		   data:	cell array containing the cell DR coordinates, labels, 
%          'movies': cell array containing movie data strack x by y by nFrames. 
%
%          Example data set here: 
%          load('/work/bioinformatics/s170480/Data/LCH/DevOps/testMovies.mat');
% Outputs:  Can manually save snapshots of plots/annotations
%
%
% Andrew R. Jamieson, Dec. 2016 updated APR 2017
%
ff = findall(0,'Tag','cellXplore');delete(ff)

ip = inputParser;
ip.KeepUnmatched = true;
ip.CaseSensitive = false;
ip.addRequired('cellDataSet', @(x) isstruct(x) || iscell(x));
ip.addOptional('movies', cell([1,length(cellDataSet)]), @iscell);
ip.addParameter('annotationSet', {}, @(x) isa(x,'containers.Map'));
ip.addParameter('DR', {}, @isstruct);
ip.parse(cellDataSet, varargin{:});
data.movies = ip.Results.movies;
data.annotationSetIn = ip.Results.annotationSet;
data.DR = ip.Results.DR; % struct {DR.PCA or DR.tSNE contains [nx2] coords} 

% Set Filter, Label, and DR Types (+ colors)
colorset = {'brgykop'};
handles.cache.plabel = {};
handles.info.zoom = false;
handles.info.GAM.state = false; 
handles.GAMfig = {};

% Initialize Label Dictionary
% check if cell array
if iscell(cellDataSet)
    disp('converning to struct array from cell array');
    cellDataSet = cell2mat(cellDataSet);
end

initializeDataStruct_Assaf();
if ~isempty(data.annotationSetIn)
    addAnnotationNoGUI(data.annotationSetIn.keys);
end

% Build out GUI/labels
initMainGUI();
initDRPanel();
initCellLabels();
handles.annoPanel={};
initFilters();
createAnnotationPanel();
createInfoPanel();
updateCustomPanels();
initCellSelect();
initMovieDisplay();
initDRAxes();

% plot everything
plotScatter;
updateCellInfo();

% Append info to UserData
cxFig = findall(0,'Tag', 'cellXplore');
cxFig.UserData.handles = handles;
cxFig.UserData.data = data;
disp('To update the UserData after events in the GUI, click the "save Notes"');
   
function initializeDataStruct_Assaf() %(MODEL)

    %     (if still a cell array use cell2mat
    [C, ia] = unique({cellDataSet.key});
    if length(C) ~= numel(cellDataSet)
        cellDataSet = cellDataSet(ia);
    end
    % master index 
    data.meta.mindex = 1:numel(cellDataSet);

    % expr date
    data.meta.expStr = {cellDataSet.expStr};
    data.meta.key = {cellDataSet.key};

    % [class labels] using struct.
    data.meta.cellType = {cellDataSet.cellType};
    data.meta.class.cellType = {cellDataSet.cellType};
    data.meta.class.metEff = {cellDataSet.metEff};
    data.meta.class.custom.all = repmat({'%+%'}, length(data.meta.mindex),1);

    % Label by experiment
    data.meta.class.MD = {cellDataSet.MD};
    data.meta.class.expr = {cellDataSet.expStr};
    
    % [notes]
    data.meta.notes = repmat({''}, length({cellDataSet.key}), 1);

    % [Annotations]
    % total set of annotation types across entire dataset.
    data.meta.anno.set = {}; %cellDataSet(1).annotations.keys; 
    data.meta.anno.tags = repmat({}, length(data.meta.mindex),1);
    % Global mapping dictionaries for entire dataset
    data.meta.anno.tagMap = containers.Map('KeyType','char','ValueType', 'any'); % tags to cells (by local master index)
    data.meta.anno.tagMapKey = containers.Map('KeyType','char','ValueType', 'any'); % tags to cells(by unique Key)
    data.meta.anno.RevTagMap = containers.Map('KeyType','char','ValueType', 'any'); % cell to tags
    % Set specific to each cell in array
    data.meta.anno.byCell = cell(1,length(data.meta.key)); % will be generated at export or save 

    % if annotations are already present for data
    if isfield(cellDataSet, 'annotations') && ...
            iscell(cellDataSet(1).annotations) && ...
            sum(cellfun(@(x) x, cellfun(@(x) ~isempty(x), {cellDataSet.annotations},'Unif',false))) >= 1
        
       % check set of annotations and update dictionaryMaps
       indx_cell = cellfun(@(x) x, cellfun(@(x) ~isempty(x), {cellDataSet.annotations},'Unif',false));
       newTags = {cellDataSet(indx_cell).annotations};
       newKeys = {cellDataSet(indx_cell).key};
%        p= containers.Map(newKeys,newTags)
       for i = 1:length(cellDataSet)
           tags = cellDataSet(i).annotations;
           key = cellDataSet(i).key;
           if ~isempty(tags)
               tagDataPointNoGUI(tags, i, key, 1);
           end
       end
    end

    % [DR types] % Ask for DR viz selection 
    feaTypes = cellDataSet(1).featureMap.keys;
    [selFeats, v] = ...
            listdlg('PromptString','Please Choose Features for DR',...
                    'SelectionMode','multiple',...
                   'ListString',feaTypes{:});
%                 feaTypes = {cellDataSet(1).featureMap.keys};
    
    %% TODO fix for multiple feat sets
    if v == 0
        selFeats = 1;
    end
    
    DRtypes = {'tSNE','PCA'};
    [selDR, v] = ...
            listdlg('PromptString','Please DR to perform',...
                    'SelectionMode','multiple',...
                    'ListString',DRtypes);
    A = cell2mat(cellfun(@(x) x(feaTypes{selFeats}), {cellDataSet.featureMap},'Unif',false))';

    if v == 0
        selDR = 2;
    end

    for dr_ = DRtypes(selDR)
        data.DR.(dr_{1}) = compute_mapping(A, dr_{1}, 2);%, 15, 30);
    end

    % [XY coord at first time point as a reference]
    data.DR.XY = [cellfun(@(x) x(1), {cellDataSet.xs}); cellfun(@(x) x(1), {cellDataSet.ys})]';
    data.extra.time = [cellfun(@(t) t(1), {cellDataSet.ts})]';

    initInfo();
end

function initInfo() % (VIEW)
    % This function translates the data to GUI handles.

    % using structs
    handles.info.labelTypes = [{'Annotate'}; fieldnames(data.meta.class)];
    handles.info.CustomTypes = fieldnames(data.meta.class.custom);

    % Set up filters
    for i=1:numel(handles.info.labelTypes)
        val = handles.info.labelTypes{i};
        if strcmp(val, 'custom') || strcmp(val, 'Annotate')  
%                 setC = unique(data.meta.class.custom);
        else
            try 
                setC = unique(data.meta.class.(val));
            catch
                setC = unique(string([data.meta.class.(val){:}]));
            end

            try 
                setC = [{'All'}, setC'];
            catch
                setC = [{'All'}, setC];
            end

            handles.info.filters.(val) = setC; 
        end

    end

    handles.info.filterTypes = fieldnames(handles.info.filters);

    % Set up filters
    for i=1:numel(handles.info.CustomTypes)
        val = handles.info.CustomTypes{i};
        setC = unique(data.meta.class.custom.(val));    
        setC = [{'All'}, setC'];
        handles.info.filters.custom.(val) = setC; 
    end

    handles.info.custfilterTypes = fieldnames(handles.info.filters.custom);

    % [DR types]
    handles.info.DRTypes_ = fieldnames(data.DR);
    % [notes]

    % [Annotations]
    handles.info.annoSet = data.meta.anno.set;
    handles.info.anno.RB = containers.Map('KeyType','char','ValueType', 'any');
    handles.selPtIdx = 1;
end

%===============================================================================
% Setup main GUI window/figure
%===============================================================================

function initMainGUI()
    xsizeF = 1375;
    ysizeF = 770;
    ysizeF = ysizeF + 5*numel(handles.info.labelTypes);
    % Create main figure
    handles.h1 = figure(...
    'Units','pixels', 'Position',[20 40 xsizeF ysizeF],...
    'Visible',get(0,'defaultfigureVisible'),...
    'Color',get(0,'defaultfigureColor'),...
    'CurrentAxesMode','manual',...
    'IntegerHandle','on',...
    'MenuBar','none',...
    'Name','cellXplore',...
    'NumberTitle','off',...
    'Tag','cellXplore',...
    'Resize','off',...
    'PaperPosition', get(0,'defaultfigurePaperPosition'),...
    'ScreenPixelsPerInchMode','manual',...
    'HandleVisibility','callback');

    handles.mainP = uipanel(...
    'Parent',handles.h1,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'Title','Cell Explorer',...
    'Tag','uipanel8',...
    'Position',[5 5 xsizeF-10 ysizeF-15],...
    'FontSize',16,...
    'FontWeight','bold');
    %-------------------------------------------------------------------------------
    % Control/Movie panels of GUI
    %-------------------------------------------------------------------------------

    % Start from bottom left
    % Left to right will be fixed.  
    % Dynamically expand vertically.
    % Standard Distance of 5 pixel between panel
    gapSize = 5;
    % Data Selection Panel
    xSizeSelectPanel = 450;
    % ---> insert dynamic size info here
    ySizeSelectPanel = 175 + 15*numel(handles.info.labelTypes); 

    % DR Panel
    xPosDR = gapSize;
    yPosDR = ySizeSelectPanel + gapSize*2;
    xSizeDRPanel= xSizeSelectPanel;
    ySizeDRPanel = 450;

    % Movie Panel
    xSizeLabelPanel = 450;

    handles.DataSel = uipanel('Parent',handles.mainP,'FontUnits','pixels','Units','pixels',...
    'Title','Data Selection/View Criterion','Tag','uipanel_select',...
    'Position',[gapSize gapSize xSizeSelectPanel ySizeSelectPanel],...
    'FontSize',13);

    handles.h2_DR = uipanel('Parent',handles.mainP, 'FontUnits','pixels', 'Units','pixels',...
    'Title','2D Visualization - Dimension Reduction',...
    'Tag','uipanel_axes',...
    'Position',[handles.DataSel.Position(1), handles.DataSel.Position(4)+handles.DataSel.Position(2)+gapSize, xSizeDRPanel, ySizeDRPanel],...
    'FontSize',13);

    handles.LabelA = uipanel('Parent',handles.mainP,'FontUnits','pixels','Units','pixels',...
    'Title','Cell Labeling',...
    'Tag','uipanel_annotate',...
    'Position',[handles.h2_DR.Position(3)+handles.h2_DR.Position(1)+gapSize, handles.DataSel.Position(2), xSizeLabelPanel 375],...
    'FontSize',13);

    handles.h_movie = uipanel(...
    'Parent',handles.mainP,'FontUnits','pixels','Units','pixels',...
    'Title','Cell Movie',...
    'Tag','uipanel_video',...
    'Position',[handles.h2_DR.Position(3)+handles.h2_DR.Position(1)+gapSize, handles.LabelA.Position(2)+handles.LabelA.Position(4)+gapSize, handles.LabelA.Position(3),...
    handles.mainP.Position(4)-handles.LabelA.Position(4)-35],...
    'FontSize',13);


    % Toggle switch for AND/OR annotation filtering


    handles.BG_AO = uibuttongroup('Parent',handles.LabelA,...
    'FontUnits','pixels','Units','pixels',...
    'Title','Filter Operation',...
    'Tag','ToggleFilter_ANDOR',...
    'Position',[325 1 100 30],...
    'FontSize', 8,...
    'SelectionChangedFcn', @updateGeneral);    

    handles.and1 = uicontrol(handles.BG_AO,'Style','radiobutton',...
                      'FontUnits','pixels','Units','pixels',...
                      'String','AND',...
                      'Position',[10 1 45 20],...
                      'HandleVisibility','off',...
                      'FontSize', 8);

    handles.or1 = uicontrol(handles.BG_AO,'Style','radiobutton',...
                      'FontUnits','pixels','Units','pixels',...
                      'String','OR',...
                      'Position',[50 1 30 20],...
                      'HandleVisibility','off','FontSize', 8);

    % Default to the OR switch    
    handles.BG_AO.SelectedObject = handles.or1;



    widthCellInfo = 330;
    handles.widthCellInfo = widthCellInfo;

    xposLabels = widthCellInfo/2;
    LabelH = 13;
    gapL = 3;
    widthString = 55;
    handles.cellInfo = uipanel(...
    'Parent',handles.mainP,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'Title','Cell Info',...
    'Tag','CellInfopanel',...
    'Position',[handles.h_movie.Position(3)+handles.h_movie.Position(1)+gapSize, handles.h_movie.Position(2),...
    widthCellInfo, numel(handles.info.labelTypes)*25+50],...
    'FontSize',13);

    % annotations & save button
    handles.notes = uicontrol(...
    'Parent',handles.LabelA,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'String','Notes here ...',...
    'Style','edit',...
    'HorizontalAlignment','left',...
    'Position',[5 5 296 22],...
    'Tag','AnnotationNotes',...
    'FontSize',13);

    % Export Full CX Data to workspace
    posE = handles.LabelA.Position;
    handles.exportData = uicontrol(...
    'Parent',handles.mainP,...
    'Style', 'pushbutton',...
    'FontUnits','pixels',...
    'Units','pixels',...
    'String','Export CellData',...
    'Position',[posE(1)+posE(3)+gapSize, posE(2), 120 25],...
    'Callback',@Export_Callback,...
    'Tag','Export CX Data',...
    'FontSize',10);
    
    function Export_Callback(varargin)
        % First save state (collect info / assign handles)
        SaveNotes_Callback(varargin{:});
        
        exportVarName = ['cellXData_' char(datetime('now','Format','ddMMMyyyy_hhmm'))];
        
        % Query for a variable name
        varname = inputdlg('Input variable name where exported cellData will be exported to',...
            'CellData export', 1, {exportVarName});
        if isempty(varname), return; end

        cxFig = findall(0,'Tag', 'cellXplore');
        assignin('base', varname{1}, cxFig.UserData.data);
        
    end
    
    % Export annotations to workspace
    posE = handles.exportData.Position;
    handles.exportDataAnno = uicontrol(...
    'Parent',handles.mainP,...
    'Style', 'pushbutton',...
    'FontUnits','pixels',...
    'Units','pixels',...
    'String','Export CellData with Annotations',...
    'Position',[posE(1)+posE(3)+gapSize, posE(2), 120 25],...
    'Callback',@ExportAnno_Callback,...
    'Tag','Export CX Data',...
    'FontSize',10);

    function ExportAnno_Callback(varargin)
        % First save state (collect info / assign handles)
        SaveNotes_Callback(varargin{:});
        
        exportVarName = ['cellMovieData_wAnno' char(datetime('now','Format','ddMMMyyyy_hhmm'))];
        
        % Query for a variable name
        varname = inputdlg('Input variable name where exported cellData will be exported to',...
            'CellData export', 1, {exportVarName});
        if isempty(varname), return; end

        % Run conversion function
        cxFig = findall(0,'Tag', 'cellXplore');
        
        cc={};
        for i=1:length(cellDataSet)
            cc{i} = cellDataSet(i);
            cc{i}.notes = {i};
            cc{i}.annotations = {};
            if isKey(data.meta.anno.RevTagMap,cc{i}.key)
                cc{i}.annotations = data.meta.anno.RevTagMap(cc{i}.key);
            end    
        end                
        
        assignin('base', varname{1}, cc);
        assignin('base', ['hashMap_' varname{1}] , data.meta.anno.tagMapKey);
        
    end
    
    
    % reconstitute original cell array with annotations
%     function re
    

    handles.SaveNotesButton = uicontrol(...
    'Parent',handles.LabelA,...
    'Style', 'pushbutton',...
    'FontUnits','pixels',...
    'Units','pixels',...
    'String','Save Notes',...
    'Position',[5 30 67 21],...
    'Callback',@SaveNotes_Callback,...
    'Tag','SaveNotes',...
    'FontSize',10);

    handles.SaveNotesNotice = uicontrol(...
    'Parent',handles.LabelA,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'Style','text',...
    'String','Saved Notes',...
    'Position',[handles.SaveNotesButton.Position(1)+handles.SaveNotesButton.Position(3)+2 handles.SaveNotesButton.Position(2) 50 13],...
    ... %     'Callback',@SaveNotes_Callback,...
    'Tag','SaveNotes',...
    'FontSize',10);
    set(handles.SaveNotesNotice, 'Visible', 'off');

    function SaveNotes_Callback(varargin)
        notes = get(handles.notes, 'String');
        data.meta.notes{handles.selPtIdx} = notes;
        set(handles.SaveNotesNotice, 'Visible', 'on');
        pause(.5);
        set(handles.SaveNotesNotice, 'Visible', 'off');

        % Append info to UserData
        cxFig = findall(0,'Tag', 'cellXplore');
        cxFig.UserData.handles = handles;
        cxFig.UserData.data = data;    
    end

    %-------------------------------------------------------------------------------
    % Movie Group Analysis
    %-------------------------------------------------------------------------------

    handles.analyzeImages = uicontrol(...
    'Parent',handles.DataSel,...
    'Style', 'pushbutton',...
    'FontUnits','pixels',...
    'Units','pixels',...
    'String','Analyze Cells',...
    'Position',[handles.DataSel.Position(3)-250 3 125 20],...
    'Callback',@analyzeGroup_Callback,...
    'Tag','SaveNotes',...
    'FontSize',12);  

    function analyzeGroup_Callback(varargin)
    % Take filter selected cells and generate new window to play
    % movies.
    %         h_ = findobj('-regexp','Tag', 'subgroupA');
    %         handles.info.GAM.fig = {}
    if isempty(handles.GAMfig)
        createAnalyzeGroupFig;
        updateAnalyzeGroupFig;
    else
        updateAnalyzeGroupFig;
    end
    end

    function createAnalyzeGroupFig()
         xsizeF = 300;
         ysizeF = 500;
        posMP = get(handles.h1, 'Position');
        % Create movieGroup
            handles.GAMfig = figure(...
            'Units','pixels',...
            'Position',[posMP(1)+10 100 xsizeF ysizeF],...
            'CurrentAxesMode','manual',...
            'IntegerHandle','on',...
            'MenuBar','none',...
            'Name','Cell Sub-Group Analysis',...
            'NumberTitle','off',...
            'Tag','subgroupA',...
            'Resize','on',...
            'CloseRequestFcn',@closeGAMfig);
    end

    function closeGAMfig(varargin)
        handles.GAMfig = {};
        handles.info.GAM.state = false;
        delete(varargin{1});
    end

    function updateAnalyzeGroupFig()
        disp('Fig already made');
        idx_f = applyFilters(handles.filters);
        handles.info.GAM.idx_f = idx_f;
        numM = length(idx_f);
        if numM > 20
            msgbox('Too Many Cells selected (must be <20)');
            close(handles.GAMfig)
            handles.GAMfig = {};
        else
            figure(handles.GAMfig);
            handles.info.GAM.state = true;
            handles.GAMs = gobjects([1 numM]);
            for i=1:numM
                if numM < 10
                    sp = subplot(numM,1,i);
                else
                    sp = subplot(ceil(numM/2),2,i);
                end
                sp.XTick = [];
                sp.YTick = [];
                sp.Box = 'off';
                sp.XColor = 'w';
                sp.YColor = 'w';
                sp.Color = [1 1 1];
                sp.Units = 'Pixels';
                handles.GAMs(i) = sp;

                imagesc(data.movies{idx_f(i)}(:,:,1), 'Parent', handles.GAMs(i), 'HitTest', 'off');
                set(handles.GAMs(i), 'XTick', []);
                set(handles.GAMs(i), 'YTick', []);
                colormap(handles.GAMs(i), gray);                     
            end
        end
        % (add legend for cells)
    end
end

%-------------------------------------------------------------------------------
% DR Type 
%-------------------------------------------------------------------------------

function initDRPanel()
    handles.DRType = uibuttongroup(...
    'Parent',handles.DataSel,...
    'FontUnits','points',...
    'Units','pixels',...
    'Title','DR View',...
    'Tag','uibuttongroup1',...
    'Position',[handles.DataSel.Position(3)-73, handles.DataSel.Position(4)-(numel(handles.info.DRTypes_)*25+25), 65, numel(handles.info.DRTypes_)*25+10],...
    'SelectionChangedFcn',@(DRType, event) DRselection(DRType, event));

    function DRselection(~, event)
       disp(['Previous: ', event.OldValue.String]);
       disp(['Current: ', event.NewValue.String]);
       disp('------------------');
       updatePlots();
    end

    handles.DRradio = gobjects(numel(handles.info.DRTypes_));
    xRB = 5;
    yRB = 8;
    for iDR=1:numel(handles.info.DRTypes_)
        handles.DRradio(iDR) = uicontrol(...
            'Parent',handles.DRType,...
            'Units','pixels',...
            'String',handles.info.DRTypes_{iDR},...
            'Style','radiobutton',...
            'Position',[xRB yRB 55 15],...
            'Tag',[handles.info.DRTypes_{iDR} '_rbutton']);
            yRB = yRB + 17;
    end

    set(handles.DRradio(1), 'Value', 1);
end

%-------------------------------------------------------------------------------
% Cell Label Menus 
%-------------------------------------------------------------------------------

function initCellLabels()
    handles.cellLabel = uicontrol(...
    'Parent',handles.LabelA,...
    'String',handles.info.labelTypes,...
    'Style','popupmenu',...
    'FontUnits','pixels',...
    'Units', 'pixels', ...
    'Value',2, ...
    'Position', [10 handles.LabelA.Position(4)-25 80 0],...
    'Callback',@updateLabel,...
    'Tag','cellLabelTypeselect',...
    'FontSize',11);

    [xCL, yCL, wCL, ~] = getPosH(handles.cellLabel);

    handles.custClass = uicontrol(...
    'Parent',handles.LabelA,...
    'String',handles.info.CustomTypes, ...
    'FontUnits','pixels',...
    'Units', 'pixels', ...
    'Style','popupmenu', ...
    'Value',1 ,...
    'Position', [xCL+wCL+5, yCL, wCL, 0],...
    'Visible', 'off',...
    'Callback',@updateLabel,...
    'Tag','customLabelMenus','FontSize',11);

    function updateLabel(source, ~)
       val = source.Value;
       maps = source.String;
       disp(['Updating Labels to : ', maps{val}]);
       disp('------------------');
       if strcmp(maps{val}, 'custom')
           set(handles.custClass, 'Visible', 'on');
       else
           set(handles.custClass, 'Visible', 'off');
       end
       updatePlots();
    end

    % Manual Label Legend
    opts = {'Parent', handles.LabelA, 'Units', 'pixels', 'Position', [11 handles.LabelA.Position(4)-150 33 83],...
            'Box' 'off','Color',[1 1 1],'XTick',[],'YTick',[]};
    axLegend = axes(opts{:});
    handles.axLegend = axLegend;
    axLegend.XColor = 'w';
    axLegend.YColor = 'w';
    % set(handles.axLegend, 'Visible', 'on', 'YAxisLocation', 'right', 'XTick', [],...
    %     'YTick', 1:8, 'YTickLabel', xlabels, 'TickLength', [0 0]);
    set(handles.axLegend, 'Visible', 'off');

    handles.dtOnOff = uicontrol(...
    'Parent',handles.LabelA,...
    'FontUnits','pixels',...
    'Units', 'pixels', ...
    'String','Show DataTips',...
    'Style','checkbox',...
    'Position',[14 137 110 17],...
    'Callback',@updateDT,...
    'Tag','checkbox1',...
    'FontSize',12, ...
    'Value', 0, ...
    'Visible', 'on');

%         handles.dtOnOff.Value = 0;

        function updateDT(source, ~)
           val = source.Value;
           disp(['Updating DataTips on/off: ', {val}]);
           disp('------------------');
           if val == 0
              set(handles.dcm_obj,'Enable','off');
           else
               dcm_obj = datacursormode(handles.h1);
               handles.dcm_obj = dcm_obj;
               set(handles.dcm_obj,'DisplayStyle','window',...
              'SnapToDataVertex','off','Enable','on');    
                      handles.dcm_obj = dcm_obj;
               set(dcm_obj,'UpdateFcn',@myupdatefcn);
           end
            updatePlots();
        end
end

%-------------------------------------------
%-------------------------------------------------------------------------------
% % Make Annotation panel
%-------------------------------------------------------------------------------

function createAnnotationPanel()
    [~, ~, w, ~] = getPosH(handles.LabelA);
    if ~isempty(handles.annoPanel)
        delete(handles.annoPanel);
    end

    numA = numel(data.meta.anno.set);

    handles.annoPanel = uipanel('Parent',handles.LabelA,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'Title','Annotations',...
    'TitlePosition', 'centertop',...
    'Tag','annoPanel',...
    'Position', [w-210, 30, 200, 55+(numA)*15], ...
    'FontSize', 12);
    xA=5; yA=5;

    handles.anno = gobjects([numA, 1]);
    handles.delAnno = gobjects([numA, 1]);
    handles.highAnno = gobjects([numA, 1]);

    % Add button for adding new annotation types
   handles.addAnno = uicontrol(...
    'Parent',handles.annoPanel,...
    'Units','pixels',...
    'FontUnits','pixels',...
    'String', 'Add New',...
    'Style', 'pushbutton',...
    'Position',[handles.annoPanel.Position(3)/2-55/2 yA 55 15],...
    'Tag','addAnnoTag',...
    'FontSize', 11,...
    'HorizontalAlignment', 'center',...
    'Callback',@addAnnotation_callback);  
    yA = yA + 15;

    for ii=1:numA
        % check tagMap
        key = data.meta.anno.set{ii};
        val = data.meta.anno.tagMap(key);
        if ismember(handles.selPtIdx,val)
           fontstyle = 'bold';
        else
           fontstyle = 'normal';
        end

        handles.anno(ii) = uicontrol(...
        'Parent',handles.annoPanel,...
        'Units','pixels',...
        'String', data.meta.anno.set{ii},...
        'Style','checkbox',...,
        'Value', ismember(handles.selPtIdx,val), ...
        'Position',[handles.annoPanel.Position(3)/2-55/2-25 yA 125 15],...
        'Tag',[data.meta.anno.set{ii} 'checkbox'],...
        'FontUnits','pixels',...
        'FontSize', 10,...
        'FontWeight', fontstyle,...
        'HorizontalAlignment', 'right',...
        'Callback', @tagDataPoint);
        if ismember(handles.selPtIdx,val)
           set(handles.anno(ii), 'BackgroundColor',[0 1 1]); 
        end

       % Add Delete button
       handles.delAnno(ii) = uicontrol(...
        'Parent',handles.annoPanel,...
        'Units','pixels',...
        'FontUnits','pixels',...
        'String', 'X',...
        'Style', 'pushbutton',...
        'Position',[handles.annoPanel.Position(3)-20 yA 15 15],...
        'Tag',[data.meta.anno.set{ii} 'delete'],...
        'FontSize', 11,...
        'HorizontalAlignment', 'center',...
        'Callback',@delAnnotation_callback);                    

        UserData = {handles.anno(ii)};
        set(handles.delAnno(ii), 'UserData', UserData);

       % Add highlight toggle
       % check previous state (use containers.Map)
       preVal = 0;
       if isKey(handles.info.anno.RB, data.meta.anno.set{ii})
            preVal = handles.info.anno.RB(data.meta.anno.set{ii});
       else
           handles.info.anno.RB(data.meta.anno.set{ii}) = 0;
       end

       handles.highAnno(ii) = uicontrol(...
        'Parent',handles.annoPanel,...
        'Units','pixels',...
        'FontUnits','pixels',...
        'String', '',...
        'Style', 'radiobutton',...
        'Value', preVal,...
        'Position',[5 yA 15 15],...
        'Tag',[data.meta.anno.set{ii} 'delete'],...
        'FontSize', 11,...
        'HorizontalAlignment', 'center',...
        'Callback',@togAnnotation_callback);                    

        if preVal == 1
            if strcmp(handles.filterAnnoMenu.String{handles.filterAnnoMenu.Value}, 'Yes')
                set(handles.highAnno(ii), 'BackgroundColor',[0 1 0]);
            end
        else
            if strcmp(handles.filterAnnoMenu.String{handles.filterAnnoMenu.Value}, 'Yes')
                set(handles.highAnno(ii), 'BackgroundColor',[1 0 0]);
            end                
        end


        UserData = {handles.anno(ii)};
        set(handles.highAnno(ii), 'UserData', UserData);
        yA = yA + 15;
    end            
        if numA >= 1
           textAnno = uicontrol(...
            'Parent',handles.annoPanel,...
            'Units','pixels',...
            'FontUnits','pixels',...
            'String', 'Tag Cell',...
            'Style', 'text',...
            'Position',[handles.annoPanel.Position(3)/2-55/2-25 yA 55 10],...
            'FontSize', 8,...
            'HorizontalAlignment', 'left');               
         textAnno = uicontrol(...
            'Parent',handles.annoPanel,...
            'Units','pixels',...
            'FontUnits','pixels',...
            'String', 'Delete',...
            'Style', 'text',...
            'Position',[handles.annoPanel.Position(3)-35 10 35 10],...
            'FontSize', 8); 

            textAnno = uicontrol(...
            'Parent',handles.annoPanel,...
            'Units','pixels',...
            'FontUnits','pixels',...
            'String', 'FilterBy',...
            'Style', 'text',...
            'Position',[1 yA 35 10],...
            'FontSize', 8,...
            'HorizontalAlignment', 'left');   
        end
end

function updateAnnotationPanel()
    numA = numel(handles.anno);

    for ih=1:numA
        hanno = handles.anno(ih);
        key = hanno.String;
        val = data.meta.anno.tagMap(key);
        if ismember(handles.selPtIdx, val)
           fontstyle = 'bold';
           set(hanno, 'BackgroundColor',[0 1 1]);    
        else
           fontstyle = 'normal';
           set(hanno, 'BackgroundColor',[.94 .94 .94]);
        end
        set(hanno, 'Value', ismember(handles.selPtIdx, val))
        set(hanno, 'FontWeight', fontstyle);
        checkRB_on(handles.highAnno(ih)); 
    end
    
end

function addAnnotation_callback(varargin)
    prompt={'Enter new Annotation Tag'};
    name = 'Add Annotation';
    default_text = {''};
    newAnnotationString = inputdlg(prompt,name,[1 30],default_text);
    if ~isempty(newAnnotationString)
        addAnnotation(newAnnotationString);
    end
end

function addAnnotation(newStr)
    if (ismember(newStr, data.meta.anno.set))
        disp('Tag already exists');
    else
        data.meta.anno.set = {data.meta.anno.set{:}, newStr{:}};
        data.meta.anno.tagMap(newStr{:}) = NaN;
        data.meta.anno.tagMapKey(newStr{:}) = {'null'};
        disp(['Added annotation tag' newStr{:}]);
        createAnnotationPanel();
    end
end

function addAnnotationNoGUI(newStrs)
    if ~iscell(newStrs)
        newStrs = {newStrs};
    end
    for newStr = newStrs
        if (ismember(newStr, data.meta.anno.set))
            disp('Tag already exists');
        else
            data.meta.anno.set = {data.meta.anno.set{:}, newStr{:}};
            data.meta.anno.tagMap(newStr{:}) = NaN;
            data.meta.anno.tagMapKey(newStr{:}) = {'null'};
            disp(['Added annotation tag [ ' newStr{:} ' ]']);
        end
    end
end


function togAnnotation_callback(src, ~)
      annoH = src.UserData{:};
      key = annoH.String;
      checkRB_on(src);
      updateFilter(handles.filterAnnoMenu, []);
end

function checkRB_on(src)
  annoH = src.UserData{:};
  key = annoH.String;
  if src.Value == src.Max
    handles.info.anno.RB(key) = 1;
    if strcmp(handles.filterAnnoMenu.String{handles.filterAnnoMenu.Value}, 'Yes')
        if strcmp(handles.BG_AO.SelectedObject.String, 'AND')
            % set as yellow (AND)
            set(src, 'BackgroundColor',[1 1 0]);    
            handles.and1.BackgroundColor = [1 1 0];
            handles.or1.BackgroundColor = [.94 .94 .94];
        else
            % set as green (OR)
            set(src, 'BackgroundColor',[0 1 0]);
            handles.BG_AO.SelectedObject.BackgroundColor = [0 1 0];
            handles.and1.BackgroundColor = [.94 .94 .94];
        end
    else
        set(src, 'BackgroundColor',[.94 .94 .94]);
    end            

  else
    handles.info.anno.RB(key) = 0;
    if strcmp(handles.filterAnnoMenu.String{handles.filterAnnoMenu.Value}, 'Yes')
        set(src, 'BackgroundColor',[1 0 0]);
    else
        set(src, 'BackgroundColor',[.94 .94 .94]);
    end

  end
end

function delAnnotation_callback(src, ~)
    annoH = src.UserData{:};
    delStr = annoH.String;
    deleteAnnotation(delStr);
end

function deleteAnnotation(tagName)
    assert(ismember(tagName, data.meta.anno.set));

    % Remove tag from cells (use KeyMap)
    rt_cells = data.meta.anno.tagMapKey(tagName);
    if length(rt_cells) > 1
    elseif length(rt_cells) == 1
        for ci = rt_cells
            % RevTagMap
            oldTags = data.meta.anno.RevTagMap(cellKey);
            newTagSet = setxor(oldTags, {tag});
            data.meta.anno.RevTagMap(cellKey) = newTagSet;
        end
        assert(strcmp(rt_cells{1}, 'null'));
    else
        error('?');
    end
    
    setOld = data.meta.anno.set;
    newSet = setxor(setOld,tagName);
    data.meta.anno.set = newSet;
    remove(data.meta.anno.tagMap, {tagName});
    remove(data.meta.anno.tagMapKey, {tagName});
    remove(handles.info.anno.RB, {tagName});
    createAnnotationPanel();
    updatePlots();
end

function tagDataPoint(src, ~)
    key = src.String;
    tag = key;
    cellKey = data.meta.key{handles.selPtIdx};
    if src.Value == src.Max % box checked
        % anno to cells
        data.meta.anno.tagMap(tag) = unique([data.meta.anno.tagMap(tag), handles.selPtIdx]);
        data.meta.anno.tagMapKey(tag) = unique([data.meta.anno.tagMapKey(tag), {cellKey}]);
        
        % cell to annos
        if isKey(data.meta.anno.RevTagMap, cellKey)
            data.meta.anno.RevTagMap(cellKey) = unique([data.meta.anno.RevTagMap(cellKey), {tag}]);
        else
            data.meta.anno.RevTagMap(cellKey) = {tag};
        end
        data.meta.anno.byCell{handles.selPtIdx} = data.meta.anno.RevTagMap(cellKey);
        
    else % box unchecked
        oldCellSet = data.meta.anno.tagMap(tag);
        newCellSet = setxor(oldCellSet, handles.selPtIdx);
        data.meta.anno.tagMap(tag) = newCellSet; 
        
        oldCellSetKey = data.meta.anno.tagMapKey(tag);
        newCellSetKey = setxor(oldCellSetKey, {cellKey});
        data.meta.anno.tagMapKey(tag) = newCellSetKey; 
        
        % RevTagMap
        oldTags = data.meta.anno.RevTagMap(cellKey);
        newTagSet = setxor(oldTags, {tag});
        data.meta.anno.RevTagMap(cellKey) = newTagSet;
        data.meta.anno.byCell{handles.selPtIdx} = data.meta.anno.RevTagMap(cellKey);

    end
    updateAnnotationPanel();
    if strcmp(handles.filterAnnoMenu.String{handles.filterAnnoMenu.Value}, 'Yes')
        updatePlots();
    end
    updateCellInfo();
end

function tagDataPointNoGUI(tags, selPtIdx, cellKey, Value)
    if ~iscell(tags) 
        tags = {tags};
    end
    for tag = tags
        tag = tag{:};
        addAnnotationNoGUI(tag); 
        if Value % box checked
            % anno to cells
            data.meta.anno.tagMap(tag) = unique([data.meta.anno.tagMap(tag), selPtIdx]);
            data.meta.anno.tagMapKey(tag) = unique([data.meta.anno.tagMapKey(tag), {cellKey}]);

            % cell to annos
            if isKey(data.meta.anno.RevTagMap, cellKey)
                data.meta.anno.RevTagMap(cellKey) = unique([data.meta.anno.RevTagMap(cellKey), {tag}]);
            else
                data.meta.anno.RevTagMap(cellKey) = {tag};
            end
            data.meta.anno.byCell{selPtIdx} = data.meta.anno.RevTagMap(cellKey);

        else % box unchecked
            oldCellSet = data.meta.anno.tagMap(tag);
            newCellSet = setxor(oldCellSet, selPtIdx);
            data.meta.anno.tagMap(tag) = newCellSet; 

            oldCellSetKey = data.meta.anno.tagMapKey(tag);
            newCellSetKey = setxor(oldCellSetKey, {cellKey});
            data.meta.anno.tagMapKey(tag) = newCellSetKey; 

            % RevTagMap
            oldTags = data.meta.anno.RevTagMap(cellKey);
            newTagSet = setxor(oldTags, {tag});
            data.meta.anno.RevTagMap(cellKey) = newTagSet;
            data.meta.anno.byCell{selPtIdx} = data.meta.anno.RevTagMap(cellKey);

        end
    end
end


%===============================================================================
% Data Label Tips / Custom Info
%===============================================================================

function createInfoPanel()
    LabelH = 13;
    widthString = 55;
    
    % Index Info Update
    ih = uicontrol(...
    'Parent',handles.cellInfo,...
    'FontUnits','pixels',...
    'HorizontalAlignment','left',...
    'Units','pixels',...
    'String','Index:',...
    'Style','text',...
    'Position',[5, 30, widthString, LabelH],...
    'Tag','text8',...
    'FontSize',10);

    handles.info.cIndex = uicontrol(...
    'Parent',handles.cellInfo,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'HorizontalAlignment','left',...
    'String',{'-'},...
    'Style','text',...
    'Position',[ih.Position(1)+ih.Position(3)+2, ih.Position(2), widthString, LabelH],...
    'Tag','IndexString',...
    'FontSize',10,...
    'FontWeight','bold');

    % Experiment Date Update
    ih = uicontrol(...
    'Parent',handles.cellInfo,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'HorizontalAlignment','left',...
    'String','ExprStr:',...
    'Style','text',...
    'Position',[ih.Position(1), ih.Position(2)+LabelH+1, widthString, 13],...
    'Tag','text14',...
    'FontSize',10);

    handles.info.ExpDateLabel_ = uicontrol(...
    'Parent',handles.cellInfo,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'HorizontalAlignment','left',...
    'String',{'01-17-2017'},...
    'Style','text',...
    'Position',[ih.Position(1)+ih.Position(3)+2, ih.Position(2), handles.widthCellInfo, LabelH],...
    'Tag','ExpDateText',...
    'FontSize',10,...
    'FontWeight','bold');

    % Cell Key ID
    ih = uicontrol(...
    'Parent',handles.cellInfo,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'HorizontalAlignment','left',...
    'String','Key:',...
    'Style','text',...
    'Position',[ih.Position(1), ih.Position(2)+LabelH+1, widthString, 13],...
    'Tag','text14',...
    'FontSize',10);

    handles.info.cellKey_ = uicontrol(...
    'Parent',handles.cellInfo,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'HorizontalAlignment','left',...
    'String',{'-'},...
    'Style','text',...
    'Position',[ih.Position(1)+ih.Position(3)+2, ih.Position(2), handles.widthCellInfo, LabelH],...
    'Tag','cellKey',...
    'FontSize',10,...
    'FontWeight','bold');

    YGapPos = ih.Position(2)+LabelH+5;

    % ==== Build out based on labels provided
    for iF=1:numel(handles.info.labelTypes)
        
        YGapPos = YGapPos+LabelH+1;

        labelName = handles.info.labelTypes{iF};

        ih = uicontrol(...
        'Parent',handles.cellInfo,...
        'FontUnits','pixels',...
        'Units','pixels',...
        'HorizontalAlignment','left',...
        'String',[labelName ': '],...
        'Style','text',...
        'Position',[ih.Position(1), YGapPos, widthString, 13],...
        'FontSize',10,...
        'Tag','text10');

        handles.info.labels.(labelName) = uicontrol(...
        'Parent',handles.cellInfo,...
        'FontUnits','pixels',...
        'Units','pixels',...
        'HorizontalAlignment','left',...
        'String',{'-'},...
        'Style','text',...
        'Position',[ih.Position(1)+ih.Position(3)+2, ih.Position(2), handles.widthCellInfo, LabelH],...
        'Tag',[labelName 'text'],...
        'FontSize',10,...
        'FontWeight','bold');      
        
    end

end

function updateCellInfo() % just a view
    
    set(handles.notes, 'String', data.meta.notes{handles.selPtIdx});

    % These will always be present
    handles.info.cIndex.String = num2str(handles.selPtIdx);
    handles.info.ExpDateLabel_.String = data.meta.expStr{handles.selPtIdx};
    handles.info.cellKey_.String = data.meta.key{handles.selPtIdx};
    handles.info.customLabel_.String = 'cust.';
    
    % These labels can change based on the class label meta data provided.
    % anything under data.meta.class (besides custom)
    for i = 1:numel(handles.info.labelTypes)
        labelName = handles.info.labelTypes{i};
        if strcmp('Annotate', labelName)
            strDisp = '';
            cellKey = data.meta.key{handles.selPtIdx};
            if isKey(data.meta.anno.RevTagMap, cellKey)
                tags = data.meta.anno.RevTagMap(cellKey);
                for tag = tags
                    if ~isempty(strDisp)
                        strDisp = [strDisp '; ' tag{:}];
                    else
                        strDisp = [strDisp tag{:}];
                    end
                end
            end
            handles.info.labels.(labelName).String = strDisp;

        elseif strcmp('custom', labelName)
        else
            handles.info.labels.(labelName).String = data.meta.class.(labelName){handles.selPtIdx};
        end
    end

end

%===============================================================================
% % ----- FIlter population 
%===============================================================================
% Filter type
% [Standard Filters]
function initFilters()
    handles.filterPanels = gobjects(numel(handles.info.filterTypes));
    yFT = 0;
    for iF=1:numel(handles.info.filterTypes)

        filterName = handles.info.filterTypes{iF};
        filterChoice = handles.info.filters.(filterName);

        handles.filterPanels(iF) = uipanel('Parent',handles.DataSel,...
         'FontUnits','pixels','Units','pixels',...
         'Title',filterName,...
         'Tag','uipanel_select',...
         'TitlePosition', 'lefttop',...
         'Position',[7 handles.DataSel.Position(4)-35-20-yFT 110 35],...
         'FontSize',10);    

        handles.filters.(filterName) = uicontrol(...
            'Parent',handles.filterPanels(iF),...
            'Units','pixels',...
            'String',filterChoice,...
            'Style','popupmenu',...
            'Position',[5 7 100 15],...
            'Callback',@updateFilter,...
            'Tag',[handles.info.filterTypes{iF} '_menuFilter'],...
            'FontWeight','bold');

        yFT = yFT + 35;

    end

    % [Filter by Annotations]
    handles.filterAnnotPanel = uipanel('Parent',handles.DataSel,...
    'FontUnits','pixels','Units','pixels',...
    'Title','Filter by Annotation',...
    'Tag','uipanel_select',...
    'TitlePosition', 'lefttop',...
    'Position',[7 handles.DataSel.Position(4)-35-20-yFT 110 35],...
    'FontSize', 10);    

    handles.filterAnnoMenu = uicontrol(...
    'Parent',handles.filterAnnotPanel,...
    'Units','pixels',...
    'String',{'No','Yes'},...
    'Style','popupmenu',...
    'Position',[5 7 100 15],...
    'Callback',@updateFilterAnno,...
    'Tag','anno_menuFilter',...
    'FontWeight','bold');


end

function updateFilterAnno(source, ~)
   val = source.Value;
   maps = source.String;
   disp(['Updating Filter by Anno to : ', maps{val}]);
   disp('------------------');
   updateAnnotationPanel();
   updatePlots();
end

% [Custom Filters]
function updateCustomPanels()

    handles.custFilterPanels = gobjects(numel(handles.info.custfilterTypes));

    if (numel(handles.info.custfilterTypes) >= 1)
        % Custom Filters 
        yFT = 0;
        for iF_ =1:numel(handles.info.custfilterTypes)

            filterName_ = handles.info.custfilterTypes{iF_};
            filterChoice_ = handles.info.filters.custom.(filterName_);

            handles.custFilterPanels(iF_) = uipanel('Parent', handles.DataSel,...
             'FontUnits','pixels','Units','pixels',...
             'Title',filterName_,...
             'Tag','uipanel_select',...
             'TitlePosition', 'lefttop',...
             'Position',[110+35 handles.DataSel.Position(4)-35-20-yFT 110 35],...
             'FontSize',10);    

            handles.custFilters.(filterName_) = uicontrol(...
                'Parent',handles.custFilterPanels(iF_),...
                'Units','pixels',...
                'String',filterChoice_,...
                'Style','popupmenu',...
                'Position',[5 7 100 15],...
                'Callback',@updateFilter,...
                'Tag',[handles.info.custfilterTypes{iF_} '_menuFilter'],...
                'FontWeight','bold');

            yFT = yFT + 35;

        end
    end

end

function updateFilter(source, ~)
   val = source.Value;
   maps = source.String;
   disp(['Updating Labels to : ', maps{val}]);
   disp('------------------');
   updatePlots();
end

% ----------------
% Manual Select Filter
% ----------------

function initCellSelect()
    handles.manualSelText = uicontrol(...
    'Parent',handles.DataSel,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'HorizontalAlignment','left',...
    'String','Select Cell Index',...
    'Style','text',...
    'Position',[handles.DataSel.Position(3)-82, 20, 80, 13],...
    'Tag','textManualIndexCellSelect',...
    'FontSize',10);

    handles.manualSel = uicontrol(...
    'Parent',handles.DataSel,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'String',arrayfun(@(x) num2str(x), 1:length(data.meta.mindex), 'UniformOutput',false), ...
    'Style','popupmenu',...
    'Value',1,...
    'Callback',@updateManSel,...
    'Position',[handles.DataSel.Position(3)-82 5 78 15],...
    'Tag','ManualIndexCellSelect',...
    'FontSize',10.667);

    handle.NextButton = uicontrol(...
    'Parent',handles.DataSel,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'HorizontalAlignment','left',...
    'String','[ Next ]',...
    'Style', 'pushbutton',...
    'Position',[handles.DataSel.Position(3)-87, 65, 40, 23],...
    'Tag','drawSelectButton',...
    'Callback',@GoToNextCell,...
    'FontSize',10);

    handle.PrevButton = uicontrol(...
    'Parent',handles.DataSel,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'BackgroundColor',[.7 .7 .7], ...
    'HorizontalAlignment','left',...
    'String','[ Prev ]',...
    'Style', 'pushbutton',...
    'Position',[handles.DataSel.Position(3)-45, 65, 40, 23],...
    'Tag','drawSelectButton',...
    'Callback',@GoToPrevCell,...
    'FontSize',10);

    handle.lassoData = uicontrol(...
    'Parent',handles.DataSel,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'HorizontalAlignment','left',...
    'String','Draw SubSet',...
    'Style', 'pushbutton',...
    'Position',[handles.DataSel.Position(3)-82, 35, 80, 23],...
    'Tag','drawSelectButton',...
    'Callback',@drawDataSelection,...
    'FontSize',10);

    function drawDataSelection(varargin)
         [selPts xc yc] = selectdata('Axes', handles.axDR);
         disp(selPts)
         % add capability to batch annotate
         % add capability to batch label?
    end

end


function GoToNextCell(varargin)
    if handles.selPtIdx == max(data.meta.mindex)
        handles.selPtIdx = 1;
    else 
        handles.selPtIdx = 1 + handles.selPtIdx;
    end
    updateAnnotationPanel();
    updatePlots();
    playMovie_GUI(); 
end

function GoToPrevCell(varargin)
    if handles.selPtIdx == 1
        handles.selPtIdx = max(data.meta.mindex);
    else
        handles.selPtIdx = -1 + handles.selPtIdx;
    end
    updateAnnotationPanel();
    updatePlots();
    playMovie_GUI();
end

             
set(handles.h1, 'KeyPressFcn', {@pb_fig, handles.h1});

function pb_fig(varargin)
    if varargin{1,2}.Character == 'n'
        GoToNextCell()
    elseif varargin{1,2}.Character == 'p'
        GoToPrevCell()
    end
end

function updateManSel(source, ~)
   val = source.Value;
   maps = source.String;
   disp(['Updating manSelect to : ', maps{val}]);
   disp(['Updating manSelect to : ', num2str(val)]);
   disp('------------------');
   handles.selPtIdx = val;
   updateAnnotationPanel();
   updatePlots();
   playMovie_GUI();
end


function updateGeneral(varargin)
   disp('Updating ...');
   disp('------------------');
   updateAnnotationPanel();
   updatePlots();
end

%===============================================================================
% Set up movie display
%===============================================================================
%-------------------------------------------------------------------------------
% Movie Panel 
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Movie Display 
%-------------------------------------------------------------------------------

function initMovieDisplay()

     h64 = uicontrol(...
    'Parent', handles.cellInfo,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'String','Play Movie',...
    'Position',[5 5 75 21],...
    'Tag','PlayMovei',...
    'FontSize',12, ...
    'Callback', @playMovie_GUI);

    h65 = uicontrol(...
    'Parent', handles.cellInfo,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'String','Master Movie',...
    'Position',[5+75+2 5 100 21],...
    'Tag','Master_Movie',...
    'FontSize',12, ...
    'Callback', @playMovieViewerAll);

    h66 = uicontrol(...
    'Parent', handles.cellInfo,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'String','Viewer [cell]',...
    'Position',[5+100+75+4 5 100 21],...
    'Tag','CellMovie',...
    'FontSize',12, ...
    'Callback', @playMovieViewerCell);
    
    opts = {'Parent', handles.h_movie, 'Units', 'pixels', 'Position',[14 36.6 handles.h_movie.Position(3)-30 handles.h_movie.Position(3)-30],...
    'Color',[1 1 1],'Box' 'off', 'XTick',[],'YTick',[]};
    axMovie = axes(opts{:});
    axMovie.XColor = 'w';
    axMovie.YColor = 'w';
    handles.axMovie = axMovie;
    handles.movies.fidx = 1; % frame index

    imagesc(data.movies{handles.selPtIdx}(:,:,1), 'Parent', handles.axMovie, 'HitTest', 'off');
    set(handles.axMovie, 'XTick', []);
    set(handles.axMovie, 'YTick', []);
    colormap(handles.axMovie, gray);

    updateSliderNF();
end

function updateSliderNF()
    % Track slider
    if ~isempty(data.movies{handles.selPtIdx})
        %  nf = size(data.movies{handles.selPtIdx},3);
        handles.movies.nf = size(data.movies{handles.selPtIdx},3);
    else
        warning('No moview provided');
        %  nf = 10; % number of frames
        handles.movies.nf = 10;
    end
    nf = handles.movies.nf;
    fidx = 1; % current frame
    handles.frameSlider = uicontrol(handles.h_movie, 'Style', 'slider', 'Units', 'pixels',...
            'Value', fidx, 'Min', 1, 'Max', nf,'SliderStep', [1/(nf-1) 0.5], ...
            'Position',[8.2 11 289.6 14],'Callback', @frameSliderRelease_Callback);   
    axMovie.Color = [1 1 1];

    addlistener(handles.frameSlider, 'Value', 'PostSet', @frameSlider_Callback);


    function frameSliderRelease_Callback(source, ~)
        val = source.Value;
        handles.movies.fidx = round(val);
        updateMovie();
        updateMovieGAM([]);
    end


    function frameSlider_Callback(~, eventdata)
        fidx_ = round(eventdata.AffectedObject.Value);
        handles.movies.fidx = round(fidx_);
        updateMovie();
        updateMovieGAM([]);
    end

end

function playMovie_GUI(varargin)
    updateSliderNF();
    handles.movies.nf = size(data.movies{handles.selPtIdx},3);
    nf = handles.movies.nf;
    tidx = handles.selPtIdx;
    i = 1;
    while (i <= nf) && (handles.selPtIdx == tidx)
        handles.movies.fidx = i;
        updateMovie();
        updateMovieGAM(i);
        pause(.05);
        i=i+1;
    end
end

function playMovieViewerCell(varargin)
end

function playMovieViewerAll(varargin)
    handles.thisMD = load(cellDataSet(handles.selPtIdx).MD);
    ff = findall(0,'Name', 'Viewer'); close(ff);
    movieViewer(handles.thisMD.MD);
    sF = findobj(0,'tag', 'slider_frame');
    sF.Value = data.extra.time(handles.selPtIdx);
    sF.Callback(sF);
    ff = findall(0,'Tag', 'viewerFig'); 
    figure(ff); 
    plot(data.DR.XY(handles.selPtIdx,1),data.DR.XY(handles.selPtIdx,2), 'or', 'markersize', 50, 'Linewidth', .25);    
end

function updateMovieGAM(Fnum)
    if ~isempty(data.movies) && handles.info.GAM.state && ~isempty(handles.GAMfig)
        idx_f = handles.info.GAM.idx_f;
        for i=1:length(idx_f)
            if ~isempty(data.movies{idx_f(i)})
                imagesc(data.movies{idx_f(i)}(:,:,handles.movies.fidx),...
                        'Parent', handles.GAMs(i), 'HitTest', 'off');
                set(handles.GAMs(i), 'XTick', []);
                set(handles.GAMs(i), 'YTick', []);
                colormap(handles.GAMs(i), gray);    
            end
        end
    end
end        

function updateMovie()
    if ~isempty(data.movies{handles.selPtIdx})
        handles.movies.nf = size(data.movies{handles.selPtIdx},3);
        imagesc(data.movies{handles.selPtIdx}(:,:,handles.movies.fidx),...
                'Parent', handles.axMovie, 'HitTest', 'off');
        set(handles.axMovie, 'XTick', []);
        set(handles.axMovie, 'YTick', []);
        colormap(handles.axMovie, gray);
    end
end        
%===============================================================================
% Set up DR viz axes
%===============================================================================

function initDRAxes()
    opts = {'Parent', handles.h2_DR, 'Units', 'pixels',...
        'Position',[15 15 handles.h2_DR.Position(3)-30 handles.h2_DR.Position(4)-65],'Color',[1 1 1],...
        'XTick',[],'YTick',[]};
    axDR = axes(opts{:});
    handles.axDR = axDR;

    % grid off;
    % Defaults
    handles.selPtIdx = 1;

    if handles.dtOnOff.Value == 1
        dcm_obj = datacursormode(handles.h1);
        handles.dcm_obj = dcm_obj;
        set(dcm_obj,'DisplayStyle','window',...
        'SnapToDataVertex','off','Enable','on');
        set(dcm_obj,'UpdateFcn',@myupdatefcn);
    end
end

%===============================================================================
% Generate Scatter Plot
%===============================================================================

function plotScatter
   
    % -----------------
    % Select Lableling
    % -----------------
        
    labeltype = handles.cellLabel.String{handles.cellLabel.Value};
    
    if strcmp(labeltype, 'custom') 
        custClass_ = handles.custClass.String{handles.custClass.Value};
        plabel = data.meta.class.custom.(custClass_);
    elseif ismember(labeltype, fieldnames(data.meta.class))
        classType_ = handles.cellLabel.String{handles.cellLabel.Value};
        plabel = data.meta.class.(classType_);
    elseif strcmp(labeltype, 'Annotate') 
        if ~isempty(handles.cache.plabel)
            plabel = handles.cache.plabel;
        else
            plabel = data.meta.class.custom.('all');    
        end
        disp('Running in annotate mode');
    else
        error('Class label not found');
    end
    
    handles.cache.plabel = plabel;

    % Generate Manual Legend
    if ~strcmp(labeltype, 'Annotate')
        try 
            [GG, GN, ~]= grp2idx(plabel);
        catch
            [GG, GN, ~]= grp2idx(cell2mat(plabel));
        end
        lcmap = cell2mat(getColors(unique(GG)));

        xlabels = cell(length(GN), 1);
        for i = 1:length(GN)
            numL = sum(ismember(GG,i));
            perctL = round(numL/length(GG),2)*100;
            xlabels{i} = [GN{i} '  (n=' num2str(numL) ') ' num2str(perctL) '%'];
        end
        
%         xlabels = GN;
        imagesc(reshape(lcmap, [size(lcmap,1) 1 3]), 'Parent', handles.axLegend);
        set(handles.axLegend, 'Visible', 'on', 'YAxisLocation', 'right', 'XTick', [],...
        'YTick', 1:8, 'YTickLabel', xlabels, 'TickLength', [0 0], 'FontSize', 7);
        set(handles.axLegend, 'Visible', 'on');
    end
    
    % get labels for plot
%     clabels = grp2idx(plabel);
    
    try 
        clabels = grp2idx(plabel);
    catch
        clabels = grp2idx(cell2mat(plabel));
    end
    
    clabels = cell2mat(getColors(clabels));
    sizeL= repmat(10,length(plabel),1);

    ji = handles.selPtIdx;
    handles.manualSel.Value = ji;  
    if handles.dtOnOff.Value == 0
        clabels(ji,:) = [0 1 1]; %[1 0 .5];
        sizeL(ji,1) = 250;
    end

    
    % ------------------------
    % Filter SubSet Data
    % ------------------------

    idx_f = applyFilters(handles.filters);
    idx_all = 1:length(data.meta.mindex);
    idx_notSel = setxor(idx_all, idx_f);
    handles.dataI = data.meta.mindex(idx_f);
    handles.dataI_ns = data.meta.mindex(idx_notSel);
    
    clabelsTemp = clabels;      
    clabelsTemp(idx_notSel, :) = repmat([1 1 1], [length(idx_notSel) 1]);
    wclabels = repmat([1 1 1], [length(idx_all) 1]);
%     end
    
    % ------------------------
    % Select DR Visualization
    % ------------------------
    
    DR_ = {handles.DRType.Children.String};
    DRtype_sel = DR_{logical([handles.DRType.Children.Value])};
    
    xyDR = data.DR.(DRtype_sel);
    X = xyDR(:,1);
    Y = xyDR(:,2);
    
    if handles.info.zoom == false

        handles.dataX = X;
        handles.dataY = Y;
        figure(handles.h1);
        handles.scat1 = scatter(handles.axDR, X, Y, sizeL-2, clabels(:,:,:),...% 'MarkerEdgeColor',clabels(:,:,:),...
            'ButtonDownFcn', @axDRCallback);
        alpha(handles.scat1, .9);
        hold on;    
        handles.scat1.SizeData(ji) = 200;
        handles.scat1.CData(ji,:) = clabels(ji,:,:);
        axis manual;

        handles.scat2 = scatter(handles.axDR, X, Y, sizeL-3, wclabels, 'filled',...
            'ButtonDownFcn', @axDRCallback);  
        handles.scat2.CData(idx_f,:) = clabels(idx_f,:,:);
        handles.scat2.CData(ji,:) = clabels(ji,:,:);
        handles.scat2.SizeData(idx_f) = 35;
        handles.scat2.SizeData(ji) = 200;
        try 
            handles.scat1.PickableParts = 'none';
            handles.scat1.HitTest = 'off';
        catch
            warning('error scat1');
        end
        
       hold off;

%        close(figure(2));    
%         end
        set(handles.axDR,'Color',[1 1 1],'Box', 'off', 'XTick',[],'YTick',[]);
        handles.axDR.Title.String = DRtype_sel;    
        handles.axDR.XColor = 'w';
        handles.axDR.YColor = 'w';
        set(handles.axDR,'Color',[1 1 1],'Box', 'off', 'XTick',[],'YTick',[]);    

    else 
        handles.dataX = X(idx_f);
        handles.dataY = Y(idx_f);
        figure(handles.h1);
        scatter(handles.axDR, X(idx_f), Y(idx_f), sizeL(idx_f), clabels(idx_f,:,:),'filled',...
            'ButtonDownFcn', @axDRCallback);
        set(handles.axDR,'Color',[1 1 1],'Box', 'off', 'XTick',[],'YTick',[]);
        handles.axDR.Title.String = DRtype_sel;    
        handles.axDR.XColor = 'w';
        handles.axDR.YColor = 'w';
        set(handles.axDR,'Color',[1 1 1],'Box', 'off', 'XTick',[],'YTick',[]);
    end
   ff = findall(0,'NumberTitle','on');delete(ff)
%    selectdata('Axes', handles.axDR);
end

%===============================================================================
% Helper functions
%===============================================================================   

function updatePlots
    updateCellInfo();
    plotScatter;
%     cxFig = findall(0,'Tag', 'cellXplore');
%     cxFig.UserData.handles = handles;
%     cxFig.UserData.data = data;
end

function txt = myupdatefcn(empt, objs)
    % handles.scat1.SizeData(handles.selPtIdx) = 10;
    idx = empt.Cursor.DataIndex;

    %% TODO - Check if filter is activated. 
%         if handles.info.zoom == false
    handles.selPtIdx = idx;
    % handles.scat1.SizeData(handles.selPtIdx) = 200;
%         else
%         if ismember(idx, handles.dataI)
% %             handles.selPtIdx = handles.dataI(idx);
%               handles.scat1.SizeData(idx) = 200;
%               handles.scat1.SizeData(idx) = 200;
%         else
% %             handles.selPtIdx = handles.dataI_ns(idx);
%         end

    %         txt = {['Index: ',num2str(handles.selPtIdx)],...
%                ['CellType: ',data.meta.cellType{handles.selPtIdx}],...
%                ['TumorType: ',data.meta.tumorTypeName{handles.selPtIdx}], ...
%                ['ExprDate :', '01-17-2017']};
    txt = {['Index: ',num2str(idx)],...
           ['xDR: ' num2str(objs.Position(1))],...
           ['yDR: ' num2str(objs.Position(2))]};

    set(handles.manualSel, 'Value', handles.selPtIdx);
    updateAnnotationPanel;
    updateCellInfo;
    playMovie_GUI();

end

function axDRCallback(varargin)
    ipt = varargin{2}.IntersectionPoint;
    a = get(gca, 'CurrentPoint');
    x0 = ipt(1,1);
    y0 = ipt(1,2);
    fx = find(round(varargin{1}.XData, 3) == round(x0,3));
    fy = find(round(varargin{1}.YData, 3) == round(y0,3));
    idx = intersect(fx,fy)
    if length(idx)>1
        %% TODO - Check if bug when inadvertently selecting multiple points.
        assert(~isempty(idx));
        idx
            idx = idx(1)
    end
%         if any(ismember(idx, handles.dataI(:)))
%             disp('Filtered');
%             handles.selPtIdx = handles.dataI(idx);
%         else
%             disp('Not Filtered');
%             handles.selPtIdx = handles.dataI_ns(idx);
%         end
    handles.selPtIdx = idx;
    updateCellInfo;
    updateAnnotationPanel;
    plotScatter; 
    playMovie_GUI;        
end

function [idx_out] = applyFilters(hinff)

    idx_out = 1:length(data.meta.mindex);
    idx_out = idx_out';
    filter_ = fieldnames(handles.filters);

    % first filter standard classes        
    for i = 1:length(filter_)

        th = hinff.(filter_{i});
        maps = th.String;  
        val = th.Value;  

        if strcmp(maps{val}, 'All')
%            disp('selecting -- all');
           idx_t = 1:length(data.meta.mindex);
           idx_t = idx_t';
        else
           idx_t = find(cellfun(@(x) strcmp(string(x), maps{val}), data.meta.class.(filter_{i}))); 
           disp(['sub-selecting ' maps{val}]);
        end
        idx_out = intersect(idx_out, idx_t);
    end
    % [then filter custom classes]
    % [then filter annotations]
    % check which annotation RadioButtons are selected
    numA = numel(data.meta.anno.set);
    if strcmp(handles.filterAnnoMenu.String{handles.filterAnnoMenu.Value}, 'Yes') && (numA >= 1)
        disp('Filtering by selected annotations');
        if strcmp(handles.BG_AO.SelectedObject.String, 'AND')
            % AND filter 
            setInx = 1:length(data.meta.mindex);
        elseif strcmp(handles.BG_AO.SelectedObject.String, 'OR')
            % OR filter 
            setInx = [];
        end
        for i=1:numA
            h_ = handles.highAnno(i);
            if (h_.Value == 1)
                annoh = h_.UserData{:};
                if strcmp(handles.BG_AO.SelectedObject.String, 'AND')
                    % AND filter 
                    setInx = intersect(setInx, data.meta.anno.tagMap(annoh.String));
                elseif strcmp(handles.BG_AO.SelectedObject.String, 'OR')
                    % OR filter
                    setInx = union(setInx, data.meta.anno.tagMap(annoh.String));
                end
            end
        end
        idx_f = unique(setInx);
        idx_out = intersect(idx_f, idx_out);
    end
end

function [RGBmat] = getColors(clabels)
   if length(unique(clabels)) <= 7
       col = colorset{:};
       RGBmat = arrayfun(@(x) let2RGB(col(x)), clabels, 'Uniform', false);
   else
       cmap = hsv(length(clabels));
       RGBmat = arrayfun(@(x) cmap(x,:), clabels, 'Uniform', false);
   end
end

function [rgbvec] = let2RGB(ltr)
    switch(lower(ltr))
        case 'r'
            rgbvec = [1 0 0];
        case 'g'
            rgbvec = [0 1 0];
        case 'b'
            rgbvec = [0 0 1];
        case 'c'
            rgbvec = [0 1 1];
        case 'm'
            rgbvec = [1 0 1];
        case 'y'
            rgbvec = [1 1 0];
        case 'w'
            rgbvec = [1 1 1];
        case 'k'
            rgbvec = [0 0 0];
        case 'o'
            rgbvec = [255/255 153/255 51/255];
        otherwise
            disp('Warning;!--colors mismatch');
    end    
end

function [x, y, w, h]=getPosH(Hin)
   x = Hin.Position(1);
   y = Hin.Position(2);
   w = Hin.Position(3); 
   h = Hin.Position(3); 
end

end




