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


function [h1] = cellXploreDR(data, varargin)


ip = inputParser;
ip.KeepUnmatched = true;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.parse(data, varargin{:});
handles.data = data;
colorset = {'brgykop'};
Celllabeltypes = { 'TumorType'; 'CellType' };
DRtypes_ = {'PCA';'tSNE'};
%===============================================================================
% Setup main GUI window/figure
%===============================================================================

% bgColor = get(0,'defaultUicontrolBackgroundColor');
% hfig = figure('Units', 'normalized', 'Position', [0.05 0.5 0.5 0.5],...
%     'PaperPositionMode', 'auto', 'Toolbar', 'figure',...
%     'Color', bgColor,...
%     'DefaultUicontrolUnits', 'pixels', 'Units', 'pixels', 'Name', ... 
%     'DRworld');
% 
% pos = get(hfig, 'Position'); % [pixels];

handles.fig11 = figure(11);
handles.fig11.NumberTitle = 'off';
handles.fig11.Name = 'Cell Movie Selector';

handles.h1 = figure(...
'Units','characters', 'Position',[135.8 45.2 97.9 39.2],...
'Visible',get(0,'defaultfigureVisible'),...
'Color',get(0,'defaultfigureColor'),...
'CurrentAxesMode','manual',...
'IntegerHandle','on',...
'MenuBar','none',...
'Name','cellXplore',...
'NumberTitle','off',...
'Tag','cellXplore',...
'Resize','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'ScreenPixelsPerInchMode','manual',...
'HandleVisibility','callback');


handles.h12 = uicontrol(...
'Parent',handles.h1,...
'FontUnits','pixels',...
'Units','characters',...
'HorizontalAlignment','left',...
'String','Cell Explorer',...
'Style','text',...
'Position',[3 36 59 2],...
'Children',[],...
'Tag','text2',...
'FontSize',18,...
'FontWeight','bold');

%-------------------------------------------------------------------------------
% Control/Movie panels of GUI
%-------------------------------------------------------------------------------

handles.h2_DR = uipanel('Parent',handles.h1, 'FontUnits','pixels', 'Units','pixels',...
'Title','2D Visualization - Dimension Reduction','Tag','uipanel_axes',...
'Position',[49 224.07 358 354],'FontSize',13,'FontSizeMode',...
get(0,'defaultuipanelFontSizeMode'));

handles.h_movie = uipanel(...
'Parent',handles.h1,'FontUnits','pixels','Units','pixels','Title','Cell Movie',...
'Tag','uipanel_video','Position',[408 377 289 202],...
'FontSize',13,'FontSizeMode',get(0,'defaultuipanelFontSizeMode'));

handles.h9 = uipanel('Parent',handles.h1,'FontUnits','pixels','Units','pixels',...
'Title','Cell Labeling','Tag','uipanel_annotate','Position',[410 77 284 293],...
'FontSize',13,'FontSizeMode',get(0,'defaultuipanelFontSizeMode'));

handles.h10 = uipanel('Parent',handles.h1,'FontUnits','pixels','Units','pixels',...
'Title','Data Selection Criterion','Tag','uipanel_select',...
'Position',[46 75 359 152],'FontSize',13,'FontSizeMode',...
get(0,'defaultuipanelFontSizeMode'));
DRtypes_

DRType = uibuttongroup(...
'Parent',handles.h10,...
'FontUnits','points',...
'Units','pixels',...
'Title','DR Type',...
'Tag','uibuttongroup1',...
'Position',[11 50 113 78],...
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

handles.DRType = DRType;

handles.h13 = uicontrol(...
'Parent',handles.DRType,...
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

handles.cellLabel = uicontrol(...
'Parent',handles.h9,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units',get(0,'defaultuicontrolUnits'),...
'String',Celllabeltypes,...
'Style','popupmenu',...
'Value',1,...
'ValueMode',get(0,'defaultuicontrolValueMode'),...
'Position',[17 243 89 22],...
'Callback',@updateLabel,...
'Tag','popupmenu1');


    function updateLabel(source, event)
       val = source.Value;
       maps = source.String;
       disp(['Updating Labels to : ', maps{val}]);
       disp('------------------');
       updatePlots()
    end

%===============================================================================
% Set up track display
%===============================================================================
nf = 10; % number of frames
fidx = 1; % current frame
handles.frameSlider = uicontrol(handles.h_movie, 'Style', 'slider', 'Units', 'pixels',...
        'Value', fidx, 'Min', 1, 'Max', nf,'SliderStep', [1/(nf-1) 0.5], ...
        'Position',[0 3 286 20]);%,'Callback', @frameSliderRelease_Callback);

%===============================================================================
% Set up DR viz axes
%===============================================================================

opts = {'Parent', handles.h2_DR, 'Units', 'pixels', 'Position', [20 20 320 290]};
axDR = axes(opts{:});
% Defaults
handles.gax = plotScatter;  datacursormode on;

% [labels lnames] = grp2idx(data.meta.cellType);
% PCA
% colorset = {'brgykop'};
% scatter(axDR, data.tSNE(:,1), data.tSNE(:,2), 14, data.meta.tumorType, 'fille');
% title_ = title('tSNE');
% set(axDR, 'Title', title_);
% set(axDR, 'Visible', 'on')
% datacursormode on;


% gax = gscatter(data.PCA(:,1), data.PCA(:,2), data.meta.cellType);



% 
% % tSNE
% title_tSNE = title('tSNE');
% scatter(haxes, x(2,:), x(4,:), 14, x(3,:));
% set(h3, 'Title', title_tSNE);




function [gax] = plotScatter
   
    ilabeltype = handles.cellLabel.Value;    
    ltyps = handles.cellLabel.String;
    labeltype = ltyps{ilabeltype};
    
    switch labeltype
        case 'CellType'
            plabel = data.meta.cellType;
        case 'TumorType'
            plabel = data.meta.tumorTypeName;
        otherwise
            plabel = data.meta.tumorTypeName;
    end

    clabels = grp2idx(plabel);
    clabels = cell2mat(getColors(clabels));

    
    iDR = find([handles.DRType.Children.Value]);
    DR_ = {handles.DRType.Children.String};
    DRtype_sel = DR_{iDR};
    
    switch DRtype_sel
       case 'PCA'
           figure(handles.h1)
           scatter(axDR, data.PCA(:,1), data.PCA(:,2), 14, clabels, 'fille');
           axDR.Title.String = 'PCA';           
           
           figure(findobj(0,'-regexp','Name', 'Movie'))
           gax = gscatter(data.PCA(:,1), data.PCA(:,2), plabel, colorset{:}, '.', 18);           
           title('PCA');
       case 'tSNE'           
           figure(handles.h1)
           scatter(axDR, data.tSNE(:,1), data.tSNE(:,2), 14, clabels, 'fille');
           axDR.Title.String = 'tSNE';
           
           figure(findobj(0,'-regexp','Name', 'Movie'))
           gax = gscatter(data.tSNE(:, 1), data.tSNE(:, 2), plabel, colorset{:}, '.', 18);
           title('tSNE');
       otherwise
            
           figure(findobj(0,'-regexp','Name', 'Movie'))
           gax = gscatter(data.PCA(:,1), data.PCA(:,2), plabel, colorset{:}, '.', 18);           
           title('PCA');
    end
end
    datacursormode on;


    function [RGBmat] = getColors(clabels)

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
       col = colorset{:}; 
       RGBmat = arrayfun(@(x) let2RGB(col(x)), clabels, 'Uniform', false);
    end

end


