function [handles] = cellXploreDR(data, varargin)
%CELLXPLORER Interactive display of movies associated with point on a 2D scatter plot.
%
% Inputs:
% 		   data:	cell array containing the cell DR coordinates, labels, 
%          'movies': cell array containing movie data strack x by y by nFrames. 
%
% Outputs:  Can manually save snapshots of plots/annotations
%
%
% Andrew R. Jamieson, Dec. 2016

ip = inputParser;
ip.KeepUnmatched = true;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParameter('movies', [], @iscell);
ip.parse(data, varargin{:});
data.movies = ip.Results.movies;

% Set Filter, Label, and DR Types (+ colors)
colorset = {'brgykop'};
handles.cache.plabel = {};
handles.info.zoom = false;
% Initialize Label Dictionary
initializeDataStruct();

    function initializeDataStruct()
        % Place commands here for standard setup based on Assaf's input
        % Standard 
        % data structure ---
        % all are indexed by mindex (respected order)
        % [DR] data.DR.{DRtype} == store 2D x,y coordinates (use fieldname)
        % [classLabels] data.meta.class.{classType} ==  store cell array of
        %   name labels for each data point
        % [custom-classLabels] data.meta.class.custom.{custClasstype} ==  store cell array of
        %   name custom labels for each data point (initialize as one class '-' for all)
        % [notes] data.meta.notes == cell array of notes for each data
        %   point
        % [annotations - types] data.meta.anno.set == cell array of annotation
        %   labels available
        % [annotations - tags] data.meta.anno.tags == array of cells contatining annotation
        %   tags (in order of mindex)
        % [annotation hash map] data.meta.anno.map == dict mapping each
        %   annotation tag as a key and storing an array of the mindex set , 
        %   used for rapid recovery of indicies and plotting. updated as
        %   new annotations are conducted.
        % 
        % Note, in annotate mode, there will be presented checkboxes and
        % editable text. -- then, there will be a button for viewing
        % highlighting the set of points that correspond to the annotated
        % selection. or a button for annotating
        % 
        % Note, in custom label mode, there will appear a new panel to the
        % left for editing the custom class categoies with menubox for selecting customclasstype, and then 
        % areas to add class labels.
        % There will always be the 'null'/'-' category, when unassigned., all
        % points are automaticall assigned this
        % note, we should try to follow the MVC mode, whereby I only access
        % data via the handles struct, and only alter data via the data
        % directly, then translate to the handle struct.
        % using containers.Map
        % dictLabels = containers.Map({'TumorType','CellType'},...
        %                          {data.meta.tumorTypeName,data.meta.cellType});
        % dictCust = containers.Map({'Custom'},...
        %                          {repmat({'none'}, length(data.meta.mindex),1)});
        % data.meta.dictLabels = dictLabels;
        % data.meta.dictCust = dictCust;
        
        % [class labels] using struct.
        data.meta.class.tumorTypeName = data.meta.tumorTypeName;
        data.meta.class.cellType = data.meta.cellType;
        data.meta.class.custom.all = repmat({'�'}, length(data.meta.mindex),1);
        % [DR types]
        data.DR.tSNE = data.tSNE;
        data.DR.PCA = data.PCA;
        % [notes]
        data.meta.notes = repmat({''}, length(data.meta.mindex),1);
        % [Annotations]
        data.meta.anno.set = {};
        data.meta.anno.tags = repmat({}, length(data.meta.mindex),1);
        data.meta.anno.tagMap = containers.Map('KeyType','char','ValueType', 'any');
        
        initInfo();
    end

    function initInfo()
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
                setC = unique(data.meta.class.(val));    
                setC = [{'All'}, setC'];
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
    end

%===============================================================================
% Setup main GUI window/figure
%===============================================================================
    function initMainGUI()
        xsizeF = 1075;
        ysizeF = 550;
        ysizeF = ysizeF + 25*numel(handles.info.labelTypes);
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
    xSizeSelectPanel = 366;
    % ---> insert dynamic size info here
    ySizeSelectPanel = 100 + 25*numel(handles.info.labelTypes); 

    % DR Panel
    xPosDR = gapSize;
    yPosDR = ySizeSelectPanel + gapSize*2;
    xSizeDRPanel= xSizeSelectPanel;
    ySizeDRPanel = 400;

    % Movie Panel
    xSizeLabelPanel = 310;

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
    'Position',[handles.h2_DR.Position(3)+handles.h2_DR.Position(1)+gapSize, handles.DataSel.Position(2), xSizeLabelPanel 330],...
    'FontSize',13);

    handles.h_movie = uipanel(...
    'Parent',handles.mainP,'FontUnits','pixels','Units','pixels',...
    'Title','Cell Movie',...
    'Tag','uipanel_video',...
    'Position',[handles.h2_DR.Position(3)+handles.h2_DR.Position(1)+gapSize, handles.LabelA.Position(2)+handles.LabelA.Position(4)+gapSize, handles.LabelA.Position(3),...
     handles.mainP.Position(4)-handles.LabelA.Position(4)-35],...
    'FontSize',13);

    widthCellInfo = 130;
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
    'Callback',@SaveNotes_Callback,...
    'Tag','SaveNotes',...
    'FontSize',10);

    set(handles.SaveNotesNotice, 'Visible', 'off');

    function SaveNotes_Callback(varargin)
        notes = get(handles.notes, 'String');
        data.meta.notes{handles.selPtIdx} = notes;
        set(handles.SaveNotesNotice, 'Visible', 'on');
        pause(.5);
        set(handles.SaveNotesNotice, 'Visible', 'off');
    end
end

initMainGUI();
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

initDRPanel();

%-------------------------------------------------------------------------------
% Cell Label Menus 
%-------------------------------------------------------------------------------
initCellLabels();

    function initCellLabels()
        handles.cellLabel = uicontrol(...
        'Parent',handles.LabelA,...
        'String',handles.info.labelTypes,...
        'Style','popupmenu',...
        'FontUnits','pixels',...
        'Units', 'pixels', ...
        'Value',1, ...
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
        'Value', 1);

            function updateDT(source, ~)
               val = source.Value;
               disp(['Updating DataTips on/off: ', {val}]);
               disp('------------------');
               if val == 0
                  set(handles.dcm_obj,'Enable','off');
               else
%                    dcm_obj = datacursormode(handles.h1);
%                    handles.dcm_obj = dcm_obj;
                   set(handles.dcm_obj,'DisplayStyle','window',...
                  'SnapToDataVertex','off','Enable','on');    
               end
                updatePlots();
            end
    end

%-------------------------------------------
%-------------------------------------------------------------------------------
% % Make Annotation panel
%-------------------------------------------------------------------------------
handles.annoPanel={};
createAnnotationPanel()

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
        'Position', [w-110, 30, 100, 55+(numA)*15], ...
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
        'Position',[xA+15 yA 55 15],...
        'Tag','addAnnoTag',...
        'FontSize', 11,...
        'HorizontalAlignment', 'center',...
        'Callback',@addAnnotation_callback);  
        yA = yA + 20;

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
            'Position',[xA+10 yA 75 15],...
            'Tag',[data.meta.anno.set{ii} 'checkbox'],...
            'FontUnits','pixels',...
            'FontSize', 10,...
            'FontWeight', fontstyle,...
            'HorizontalAlignment', 'right',...
            'Callback', @tagDataPoint);

           % Add Delete button
           handles.delAnno(ii) = uicontrol(...
            'Parent',handles.annoPanel,...
            'Units','pixels',...
            'FontUnits','pixels',...
            'String', 'X',...
            'Style', 'pushbutton',...
            'Position',[xA+75 yA 15 15],...
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
            'Position',[xA-5 yA 15 15],...
            'Tag',[data.meta.anno.set{ii} 'delete'],...
            'FontSize', 11,...
            'HorizontalAlignment', 'center',...
            'Callback',@togAnnotation_callback);                    
              
            UserData = {handles.anno(ii)};
            set(handles.highAnno(ii), 'UserData', UserData);
            yA = yA + 15;
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
            else
               fontstyle = 'normal';
            end
            set(hanno, 'Value', ismember(handles.selPtIdx, val))
            set(hanno, 'FontWeight', fontstyle);    
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
            data.meta.anno.tagMap(newStr{:}) = 0;
            disp(['Added annotation tag' newStr{:}]);
            createAnnotationPanel();
        end
    end

    function togAnnotation_callback(src, ~)
          annoH = src.UserData{:};
          key = annoH.String;
          if src.Value == src.Max
            handles.info.anno.RB(key) = 1;
          else
            handles.info.anno.RB(key) = 0;  
          end
          updateFilter(handles.filterAnnoMenu, []);
    end

    function delAnnotation_callback(src, ~)
        annoH = src.UserData{:};
        delStr = annoH.String;
        deleteAnnotation(delStr);
    end

    function deleteAnnotation(tagName)
        assert(ismember(tagName, data.meta.anno.set));
        setOld = data.meta.anno.set;
        newSet = setxor(setOld,tagName);
        data.meta.anno.set = newSet;
        remove(data.meta.anno.tagMap, {tagName});
        remove(handles.info.anno.RB, {tagName});
        createAnnotationPanel();
    end

    function tagDataPoint(src, ~)
        key = src.String;
        if src.Value == src.Max
            data.meta.anno.tagMap(key) = unique([data.meta.anno.tagMap(key), handles.selPtIdx]);
%             set(src, 'FontWeight', 'Bold');
        else
            oldSet = data.meta.anno.tagMap(key);
            newSet = setxor(oldSet, handles.selPtIdx);
            data.meta.anno.tagMap(key) = newSet;    
%             set(src, 'FontWeight', 'Normal');
        end
        updateAnnotationPanel();
        if strcmp(handles.filterAnnoMenu.String{handles.filterAnnoMenu.Value}, 'Yes')
            updatePlots();
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

    % TumorType Update
    ih = uicontrol(...
    'Parent',handles.cellInfo,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'HorizontalAlignment','left',...
    'String','TumorType:',...
    'Style','text',...
    'Position',[ih.Position(1), ih.Position(2)+LabelH+1, widthString, 13],...
    'FontSize',10,...
    'Tag','text10');

    handles.info.TumorTypeLabel_ = uicontrol(...
    'Parent',handles.cellInfo,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'HorizontalAlignment','left',...
    'String',{'-T'},...
    'Style','text',...
    'Position',[ih.Position(1)+ih.Position(3)+2, ih.Position(2), widthString, LabelH],...
    'Tag','TumorTypeText',...
    'FontSize',10,...
    'FontWeight','bold');

    % CellType Update
    ih = uicontrol(...
    'Parent',handles.cellInfo,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'HorizontalAlignment','left',...
    'String','CellType:',...
    'Style','text',...
    'Position',[ih.Position(1), ih.Position(2)+LabelH+1, widthString, 13],...
    'Tag','text12',...
    'FontSize',10);

    handles.info.CellTypeLabel_ = uicontrol(...
    'Parent',handles.cellInfo,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'HorizontalAlignment','left',...
    'String',{'-'},...
    'Style','text',...
    'Position',[ih.Position(1)+ih.Position(3)+2, ih.Position(2), widthString, LabelH],...
    'Tag','CellTypeText',...
    'FontSize',10,...
    'FontWeight','bold');

    % Experiment Date Update
    ih = uicontrol(...
    'Parent',handles.cellInfo,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'HorizontalAlignment','left',...
    'String','ExprDate:',...
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
    'Position',[ih.Position(1)+ih.Position(3)+2, ih.Position(2), widthString, LabelH],...
    'Tag','ExpDateText',...
    'FontSize',10,...
    'FontWeight','bold');
    % 
    % h62 = uicontrol(...
    % 'Parent',handles.cellInfo,...
    % 'FontUnits','pixels',...
    % 'Units','pixels',...
    % 'Position',[14 54 63 13],...
    % 'HorizontalAlignment','right',...
    % 'String','Custom:',...
    % 'Style','text',...
    % 'Tag','text16',...
    % 'FontSize',10);
    % 
    % handles.info.customLabel_ = uicontrol(...
    % 'Parent',handles.cellInfo,...
    % 'FontUnits','pixels',...
    % 'Units','pixel',...
    % 'HorizontalAlignment','left',...
    % 'String',{'-'},...
    % 'Style','text',...
    % 'Position',[xposLabels+5 54 80 13],...
    % 'Tag','CustomText',...
    % 'FontWeight','bold',...
    % 'FontSize',10);

    h64 = uicontrol(...
    'Parent', handles.cellInfo,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'String','Play Movie',...
    'Position',[5 5 75 21],...
    'Tag','PlayMovei',...
    'FontSize',12, ...
    'Callback', @playMovie_GUI);
end
createInfoPanel();

%===============================================================================
% % ----- FIlter population 
%===============================================================================

% Filter type
% [Standard Filters]

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

function updateFilterAnno(source, ~)
   val = source.Value;
   maps = source.String;
   disp(['Updating Filter by Anno to : ', maps{val}]);
   disp('------------------');
   updatePlots();
end

updateCustomPanels();

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
                 'Position',[110+12 handles.DataSel.Position(4)-35-20-yFT 110 35],...
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

initCellSelect();

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


%===============================================================================
% Set up movie display
%===============================================================================

%-------------------------------------------------------------------------------
% Movie Panel 
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Movie Display 
%-------------------------------------------------------------------------------
% initialize Movie
initMovieDisplay();
    
    function initMovieDisplay()

        opts = {'Parent', handles.h_movie, 'Units', 'pixels', 'Position',[14 36.6 handles.h_movie.Position(3)-30 206.8],...
        'Color',[1 1 1],'Box' 'off', 'XTick',[],'YTick',[]};
        axMovie = axes(opts{:});
        axMovie.XColor = 'w';
        axMovie.YColor = 'w';
        handles.axMovie = axMovie;
        handles.movies.fidx = 1; % frame index
        
        imagesc(data.movies{1}(:,:,1), 'Parent', handles.axMovie, 'HitTest', 'off');
        set(handles.axMovie, 'XTick', []);
        set(handles.axMovie, 'YTick', []);
        colormap(handles.axMovie, gray);


        % Track slider
        if ~isempty(data.movies)
            nf = size(data.movies{1},3);
        else
            warning('No moview provided');
            nf = 10; % number of frames
        end

        handles.movies.nf = nf;

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
        end

    
        function frameSlider_Callback(~, eventdata)
            fidx_ = round(eventdata.AffectedObject.Value);
            handles.movies.fidx = round(fidx_);
            updateMovie();
        end
    end

    function playMovie_GUI(varargin)
        nf = handles.movies.nf;
        tidx = handles.selPtIdx;
        i = 1;
        while (i <= nf) && (handles.selPtIdx == tidx)
            handles.movies.fidx = i;
            updateMovie();
            pause(.05);
            i=i+1;
        end
    end

    function updateMovie()
        if ~isempty(data.movies)
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

initDRAxes();

    function initDRAxes()
        opts = {'Parent', handles.h2_DR, 'Units', 'pixels',...
            'Position',[15 15 handles.h2_DR.Position(3)-30 handles.h2_DR.Position(4)-65],'Color',[1 1 1],...
            'XTick',[],'YTick',[]};
        axDR = axes(opts{:});
        handles.axDR = axDR;

        % grid off;
        % Defaults
        handles.selPtIdx = 1;

        dcm_obj = datacursormode(handles.h1);
        handles.dcm_obj = dcm_obj;
        set(dcm_obj,'DisplayStyle','window',...
        'SnapToDataVertex','off','Enable','on');
        set(dcm_obj,'UpdateFcn',@myupdatefcn);
    end

% plot everything
plotScatter;

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
        [GG, GN, ~]= grp2idx(plabel);
        lcmap = cell2mat(getColors(unique(GG)));
        xlabels = GN;
        imagesc(reshape(lcmap, [size(lcmap,1) 1 3]), 'Parent', handles.axLegend);
        set(handles.axLegend, 'Visible', 'on', 'YAxisLocation', 'right', 'XTick', [],...
        'YTick', 1:8, 'YTickLabel', xlabels, 'TickLength', [0 0]);
        set(handles.axLegend, 'Visible', 'on');
    end
    
    % get labels for plot
    clabels = grp2idx(plabel);
    clabels = cell2mat(getColors(clabels));
    sizeL= repmat(12,length(plabel),1);

    ji = handles.selPtIdx;
    handles.manualSel.Value = ji;  
    if handles.dtOnOff.Value == 0
        clabels(ji,:) = [0 1 1]; %[1 0 .5];
        sizeL(ji,1) = 200;
    end

    % ------------------------
    % Filter SubSet Data
    % ------------------------

    idx_f = applyFilters(handles.filters);
    handles.dataI = data.meta.mindex(idx_f);
    
%     % Color non-selected points white
%     idx_all = 1:length(data.meta.mindex);
%     idx_notSel = setxor(idx_all, idx_f);
%     if ~isempty(idx_notSel)
%         clabels(idx_notSel, :) = repmat([1 1 1], [length(idx_notSel) 1]);
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
        s = scatter(handles.axDR, X, Y, sizeL, clabels(:,:,:),'filled',...
            'ButtonDownFcn', @axDRCallback);
        if length(idx_f) ~= length(data.meta.mindex) 
            alpha(s, .20);
            hold on;
            axis manual;
            scatter(handles.axDR, X(idx_f), Y(idx_f), 25, clabels(idx_f,:,:),'filled',...
                'ButtonDownFcn', @axDRCallback);
            hold off;
        end
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
end

%===============================================================================
% Helper functions
%===============================================================================   
    function updateCellInfo()
        set(handles.notes, 'String', data.meta.notes{handles.selPtIdx});

        % ------------
        % update Cell Info Panel
        handles.info.cIndex.String = num2str(handles.selPtIdx);
        handles.info.TumorTypeLabel_.String = data.meta.tumorTypeName{handles.selPtIdx};
        handles.info.CellTypeLabel_.String = data.meta.cellType{handles.selPtIdx}; 
        handles.info.ExpDateLabel_.String = '11-11-2017';
        handles.info.customLabel_.String = 'cust.';
        
        % ------------        
    end
        
    function updatePlots
        updateCellInfo();
%         updateAnnotationPanel();
        plotScatter;
    end
    
    function txt = myupdatefcn(empt, objs)
        % Customizes text of data tips
        idx = empt.Cursor.DataIndex;
        if handles.info.zoom == false
            handles.selPtIdx = idx;
        else
            handles.selPtIdx = handles.dataI(idx);
        end
%         txt = {['Index: ',num2str(handles.selPtIdx)],...
%                ['CellType: ',data.meta.cellType{handles.selPtIdx}],...
%                ['TumorType: ',data.meta.tumorTypeName{handles.selPtIdx}], ...
%                ['ExprDate :', '01-17-2017']};
        txt = {['xDR: ' num2str(objs.Position(1))],...
               ['yDR: ' num2str(objs.Position(2))]};

        set(handles.manualSel, 'Value', handles.selPtIdx);
        updateAnnotationPanel;
        updateCellInfo;
        playMovie_GUI();
    end

    function axDRCallback(varargin)
        ipt = varargin{2}.IntersectionPoint;
        x0 = ipt(1,1);
        y0 = ipt(1,2);
        fx = find(round(varargin{1}.XData, 5) == round(x0,5));
        fy = find(round(varargin{1}.YData, 5) == round(y0,5));
        idx = intersect(fx,fy);
        if handles.info.zoom == false
            handles.selPtIdx = idx;
        else
            handles.selPtIdx = handles.dataI(idx);
        end
        if length(handles.selPtIdx) > 1
            handles.selPtIdx = handles.selPtIdx(1);
        end
     
        plotScatter; 
        updateCellInfo;
        updateAnnotationPanel;
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
               disp('selecting -- all');
               idx_t = 1:length(data.meta.mindex);
               idx_t = idx_t';
            else
               idx_t = find(cellfun(@(x) strcmp(x, maps{val}), data.meta.class.(filter_{i}))); 
               disp(['sub-selecting ' maps{val}]);
            end
            idx_out = intersect(idx_out, idx_t);
        end
        
        % then filter custom classes
        
        % then filter annotations
        % check which annotation RadioButtons are selected
        numA = numel(data.meta.anno.set);
        if strcmp(handles.filterAnnoMenu.String{handles.filterAnnoMenu.Value}, 'Yes') && (numA >= 1)
            disp('Filtering by selected annotations');
            setInx = 1:length(data.meta.mindex);
            for i=1:numA
                h_ = handles.highAnno(i);
                if (h_.Value == 1)
                    annoh = h_.UserData{:};
                    setInx = intersect(setInx, data.meta.anno.tagMap(annoh.String));
%                     setInx = [data.meta.anno.tagMap(annoh.String), setInx]; 
                end
            end
            idx_f = unique(setInx);
            idx_out = intersect(idx_f, idx_out);
        end
        
    end

    function [RGBmat] = getColors(clabels)
       col = colorset{:}; 
       RGBmat = arrayfun(@(x) let2RGB(col(x)), clabels, 'Uniform', false);
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




