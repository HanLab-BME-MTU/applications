%cellXploreMan(data, varargin) interactive display of movies via 2D DR plot.
%
% Inputs:
% 		  data:	cell array containing the cell DR coordinates, labels, 
%			% and movie paths as well as other metadata.
%
% Outputs:  		Can manually save snapshots of plots/annotations
%
%
% Andrew R. Jamieson, Dec. 2016


function [handles] = cellXploreDR(data, varargin)


ip = inputParser;
ip.KeepUnmatched = true;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.parse(data, varargin{:});
handles.data = data;

% Set Filter, Label, and DR Types (+ colors)
colorset = {'brgykop'};
labelTypes = { 'TumorType'; 'CellType' };
% TumorTypeLabels = {'Malignant'; 'Benign'; 'All'};
[G, G2] = grp2idx(data.meta.cellType);
[Gi, G2i] = grp2idx(data.meta.tumorTypeName);
cellTypes = [{ 'All' }, G2']; 
TumorTypeLabels = [{ 'All' }, G2i']; 
DRtypes_ = {'PCA';'tSNE'};

handles.info.DRtypes_ = DRtypes_;
handles.info.cellTypes = cellTypes;
handles.info.TumorTypeLabels = TumorTypeLabels;

%===============================================================================
% Setup main GUI window/figure
%===============================================================================

% Separate gscatter plot for cursor
% handles.fig11 = figure(11);
% handles.fig11.NumberTitle = 'off';
% handles.fig11.Name = 'Cell Movie Selector';

% Create main figure
handles.h1 = figure(...
'Units','pixels', 'Position',[10 10 735 672],...
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

% Main Title
handles.h12 = uicontrol(...
'Parent',handles.h1,...
'FontUnits','pixels',...
'Units','pixels',...
'HorizontalAlignment','left',...
'String','Cell Explorer',...
'Style','text',...
'Position',[20.6 643.4 306.8 23.2],...
'Children',[],...
'Tag','text2',...
'FontSize',18,...
'FontWeight','bold');

%-------------------------------------------------------------------------------
% Control/Movie panels of GUI
%-------------------------------------------------------------------------------

handles.h2_DR = uipanel('Parent',handles.h1, 'FontUnits','pixels', 'Units','pixels',...
'Title','2D Visualization - Dimension Reduction','Tag','uipanel_axes',...
'Position',[14.6 237.8 366.4 401.2],'FontSize',13,'FontSizeMode',...
get(0,'defaultuipanelFontSizeMode'));

handles.h_movie = uipanel(...
'Parent',handles.h1,'FontUnits','pixels','Units','pixels','Title','Cell Movie',...
'Tag','uipanel_video','Position',[387.4 376.6 308 262.8],...
'FontSize',13,'FontSizeMode',get(0,'defaultuipanelFontSizeMode'));

handles.LabelA = uipanel('Parent',handles.h1,'FontUnits','pixels','Units','pixels',...
'Title','Cell Labeling','Tag','uipanel_annotate','Position',[386.6 76.6 307.6 299.6],...
'FontSize',13,'FontSizeMode',get(0,'defaultuipanelFontSizeMode'));

handles.DataSel = uipanel('Parent',handles.h1,'FontUnits','pixels','Units','pixels',...
'Title','Data Selection Criterion','Tag','uipanel_select',...
'Position',[15.4 77.4 365.2 159.2],'FontSize',13,'FontSizeMode',...
get(0,'defaultuipanelFontSizeMode'));

DRType = uibuttongroup(...
'Parent',handles.DataSel,...
'FontUnits','points',...
'Units','pixels',...
'Title','DR Type',...
'Tag','uibuttongroup1',...
'Position',[240.6 60.4 113.2 78],...
'SelectionChangedFcn',@(DRType, event) DRselection(DRType, event));

function DRselection(bg, event)
   disp(['Previous: ', event.OldValue.String]);
   disp(['Current: ', event.NewValue.String]);
   disp('------------------');
   updatePlots()
%        setappdata
end

function updatePlots
    plotScatter;
end
    
%-------------------------------------------------------------------------------
% DR Type 
%-------------------------------------------------------------------------------
handles.DRType = DRType;

handles.h13 = uicontrol(...
'Parent',handles.DRType,...cellLabel
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','pixels',...
'String',DRtypes_{2},...
'Style','radiobutton',...
'Value',1,...
'Position',[11 35 80 17],...
'Tag','radiobutton1');

handles.h14 = uicontrol(...
'Parent',handles.DRType,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','pixels',...
'String',DRtypes_{1},...
'Style','radiobutton',...
'Position',[12 10 80 17],...
'Tag','radiobutton2');

%-------------------------------------------------------------------------------
% Cell Label Menus 
%-------------------------------------------------------------------------------

handles.cellLabel = uicontrol(...
'Parent',handles.LabelA,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units',get(0,'defaultuicontrolUnits'),...
'String',labelTypes,...
'Style','popupmenu',...
'Value',1,...
'ValueMode',get(0,'defaultuicontrolValueMode'),...
'Position',[10 257 89 22],...
'Callback',@updateLabel,...
'Tag','cellLabelTypeselect');

function updateLabel(source, event)
   val = source.Value;
   maps = source.String;
   disp(['Updating Labels to : ', maps{val}]);
   disp('------------------');
   updatePlots()
end


% Manual Label Legend
opts = {'Parent', handles.LabelA, 'Units', 'pixels', 'Position', [11 161 33 83],...
        'Box' 'off','Color',[1 1 1],'XTick',[],'YTick',[]};
axLegend = axes(opts{:});
handles.axLegend = axLegend;
axLegend.XColor = 'w';
axLegend.YColor = 'w';
% set(handles.axLegend, 'Visible', 'on', 'YAxisLocation', 'right', 'XTick', [],...
%     'YTick', 1:8, 'YTickLabel', xlabels, 'TickLength', [0 0]);
set(handles.axLegend, 'Visible', 'off');
%===============================================================================
% % ----- FIlter population 
%===============================================================================

% Main Sub-menu title
handles.filterText = uicontrol(...
'Parent',handles.DataSel,...
'FontUnits','pixels',...
'Units','pixels',...
'HorizontalAlignment','left',...
'String','Population Subset',...
'Style','text',...
'Position',[9 120.6 109.2 17.2],...
'Tag','text3',...
'FontSize',13);

% Filter type
handles.filterTextT = uicontrol(...
'Parent',handles.DataSel,...
'FontUnits','pixels',...
'String','TumorType',...
'Style','text',...
'Position',[99.8 101 60.8 13.2],...
'Tag','text4',...
'FontSize',10.66);

handles.filters.tumorTypeName = uicontrol(...
'Parent',handles.DataSel,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units',get(0,'defaultuicontrolUnits'),...
'String',TumorTypeLabels,...
'Style','popupmenu',...
'Value',1, ...
'ValueMode',get(0,'defaultuicontrolValueMode'),...
'Position',[10.6 96.6 89.2 22],...
'Callback',@updateFilter,...
'Tag','popupmenu2');

    function updateFilter(source, event)
       val = source.Value;
       maps = source.String;
       disp(['Updating Labels to : ', maps{val}]);
       disp('------------------');
       updatePlots()
    end


% CellType Filter

handles.filtersTextC = uicontrol(...
'Parent',handles.DataSel,...
'FontUnits','pixels',...
'Units','pixels',...
'String','CellType',...
'Style','text',...
'Position',[100 80.6 53.6 13.2],...
'Tag','CellTypeText',...
'FontSize',10.66);

handles.filters.cellType = uicontrol(...
'Parent',handles.DataSel,...
'FontUnits','pixels',...
'Units','pixels',...
'String',cellTypes,...
'Style','popupmenu',...
'Value',1, ...
'Position',[10.6 75.8 89.2 22],...
'Callback',@updateFilterC,...
'Tag','popupmenu2');

function updateFilterC(source, event)
   val = source.Value;
   maps = source.String;
   disp(['Updating Labels to : ', maps{val}]);
   disp('------------------');
   updatePlots()
end

% ----------------
% CellType Filter
% ----------------

h21 = uicontrol(...
'Parent',handles.DataSel,...
'FontUnits','pixels',...
'Units','pixels',...
'HorizontalAlignment','left',...
'String','Select Cell Index',...
'Style','text',...
'Position',[10.6 30.2 107.2 13.2],...
'Tag','textManualIndexCellSelect',...
'FontSize',10.6);

h22 = uicontrol(...
'Parent',handles.DataSel,...
'FontUnits','pixels',...
'Units','pixels',...
'String',{  '1','2','3' },...
'Style','popupmenu',...
'Value',1,...
'Position',[10.6 10.6 78 15.6],...
'Tag','ManualIndexCellSelect',...
'FontSize',10.667);


%===============================================================================
% Set up movie display
%===============================================================================
opts = {'Parent', handles.h_movie, 'Units', 'pixels', 'Position',[18.2 36.6 275.6 206.8],...
    'Color',[1 1 1],'Box' 'off', 'XTick',[],'YTick',[]};
axMovie = axes(opts{:});
axMovie.XColor = 'w';
axMovie.YColor = 'w';
% Track slider
nf = 10; % number of frames
fidx = 1; % current frame
handles.frameSlider = uicontrol(handles.h_movie, 'Style', 'slider', 'Units', 'pixels',...
        'Value', fidx, 'Min', 1, 'Max', nf,'SliderStep', [1/(nf-1) 0.5], ...
        'Position',[8.2 11 289.6 14]);%,'Callback', @frameSliderRelease_Callback);   
axMovie.Color = [1 1 1];


%===============================================================================
% Set up DR viz axes
%===============================================================================

opts = {'Parent', handles.h2_DR, 'Units', 'pixels', 'Position',[20 40.2 323.2 308],'Color',[1 1 1],...
    'XTick',[],'YTick',[]};
axDR = axes(opts{:});
handles.axDR = axDR;

% grid off;
% Defaults
handles.selPtIdx = 1;
dcm_obj = datacursormode(handles.h1);
set(dcm_obj,'DisplayStyle','datatip',...
'SnapToDataVertex','off','Enable','on');
set(dcm_obj,'UpdateFcn',@myupdatefcn)
% c_info = getCursorInfo(dcm_obj);
    
handles.gax = plotScatter;  
% datacursormode on;


%===============================================================================
% Generate Scatter Plot
%===============================================================================

function [gax] = plotScatter
   
    ilabeltype = handles.cellLabel.Value;    
    ltyps = handles.cellLabel.String;
    labeltype = ltyps{ilabeltype};
    gax = [];

    % -----------------
    % Select Lableling
    % -----------------

    switch labeltype
        case 'CellType'
            plabel = data.meta.cellType;
        case 'TumorType'
            plabel = data.meta.tumorTypeName;
        otherwise
            plabel = data.meta.tumorTypeName;
    end
    
    ji = handles.selPtIdx;
    
    % Generate Manual Legend
    [GG GN  GL]= grp2idx(plabel);
    lcmap = cell2mat(getColors(unique(GG)));
    xlabels = GN;
    imagesc(reshape(lcmap, [size(lcmap,1) 1 3]), 'Parent', handles.axLegend);
    set(handles.axLegend, 'Visible', 'on', 'YAxisLocation', 'right', 'XTick', [],...
    'YTick', 1:8, 'YTickLabel', xlabels, 'TickLength', [0 0]);
    set(handles.axLegend, 'Visible', 'on');

    % get labels for plot
    clabels = grp2idx(plabel);
    clabels = cell2mat(getColors(clabels));
    clabels(ji,:) = [0 1 1];%[1 0 .5];
    
        
    sizeL= repmat(12,length(plabel),1);
    sizeL(ji,1) = 100;
    
    % ------------------------
    % Filter SubSet Data
    % ------------------------

    idx_f = applyFilters(handles.filters);
    
    % ------------------------
    % Select DR Visualization
    % ------------------------
    
    iDR = find([handles.DRType.Children.Value]);
    DR_ = {handles.DRType.Children.String};
    DRtype_sel = DR_{iDR};
    
    switch DRtype_sel
       case 'PCA'
           figure(handles.h1);
           scatter(axDR, data.PCA(idx_f,1), data.PCA(idx_f,2), sizeL(idx_f), clabels(idx_f,:,:),'filled');%,'ButtonDownFcn', @axDRCallback);
           set(axDR,'Color',[1 1 1],'Box', 'off', 'XTick',[],'YTick',[]);
           axDR.Title.String = 'PCA';           
       case 'tSNE'           
           figure(handles.h1);
           scatter(axDR, data.tSNE(idx_f,1), data.tSNE(idx_f,2), sizeL(idx_f), clabels(idx_f,:,:), 'filled');%, 'ButtonDownFcn', @axDRCallback);
           axDR.Title.String = 'tSNE';
        otherwise
    end
    axDR.XColor = 'w';
    axDR.YColor = 'w';
    set(axDR,'Color',[1 1 1],'Box', 'off', 'XTick',[],'YTick',[]);
end


    
    
%===============================================================================
% Helper functions
%===============================================================================   
%  

    function txt = myupdatefcn(empt, event_obj)
        % Customizes text of data tips
        handles.selPtIdx = empt.Cursor.DataIndex;
        pos = get(event_obj,'Position');
        txt = {['Time: ',num2str(pos(1))],...
                  ['Amplitude: ',num2str(pos(2))]};
      plotScatter;
    end

    function axDRCallback(varargin)
        a = get(gca, 'CurrentPoint');
        x0 = a(1,1);
        y0 = a(1,2);
        x0
        y0
      iDR = find([handles.DRType.Children.Value]);
      DR_ = {handles.DRType.Children.String};
      DRtype_sel = DR_{iDR};
      switch DRtype_sel
       case 'PCA'
       case 'tSNE'           
       otherwise
      end
        
    end

    function [idx_out] = applyFilters(hinff)

        fc = fieldnames(hinff);
        idx_out = 1:length(data.meta.mindex);
        idx_out = idx_out';

        for i = 1:length(fc)

            th = hinff.(fc{i});
            maps = th.String;  
            val = th.Value;  

            if strcmp(maps{val}, 'All')
               disp('selecting -- all');
               idx_t = 1:length(data.meta.mindex);
               idx_t = idx_t';
            else
               idx_t = find(cellfun(@(x) strcmp(x, maps{val}), data.meta.(fc{i}))); 
               disp(['sub-selecting ' maps{val}]);
            end
            idx_out = intersect(idx_out, idx_t);
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
end


