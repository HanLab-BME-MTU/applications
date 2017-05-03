function [cxFig handles] = dopeAnnotator(cellDataSet, varargin)
%DOPEANNOTATOR Interactive display of movies associated with point on a 2D scatter plot.
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
ff = findall(0,'Tag','dopeAnnotator');delete(ff)
data.movies = {};
data.annotationSetIn = {};

if nargin >= 1
    ip = inputParser;
    ip.KeepUnmatched = true;
    ip.CaseSensitive = false;
    ip.addRequired('cellDataSet', @(x) isstruct(x) || iscell(x) || exist(x, 'file')==2 || @(x) isa(x, 'MovieList'));
    ip.addParameter('movies', cell([1,length(cellDataSet)]), @iscell);
    ip.addParameter('annotationSet', {}, @(x) isa(x,'containers.Map'));
    ip.addParameter('DR', {}, @isstruct);
    ip.parse(cellDataSet, varargin{:});
    data.movies = ip.Results.movies;
    data.annotationSetIn = ip.Results.annotationSet;
    data.DR = ip.Results.DR; % struct {DR.PCA or DR.tSNE contains [nx2] coords} 
else
    cellDataSet = {};
end

data.moviesPath = 'none provided...';
handles.FastPlotMode = false;

handles.cache.plabel = {};
handles.info.zoom = false;
handles.info.GAM.state = false; 
handles.GAMfig = {};
handles.pointSize = 4;
handles.pointSizeFilter = 20;
handles.zoomDR = 'off';
handles.shadowPoints = true;
handles.forceUpdate = false;
handles.choiceLoadMD = 'no';
handles.flyLoad = false; % load MovieDatas on the fly
handles.logfile = '/work/bioinformatics/shared/dope/export/.AnnotationsLog.txt';
handles.backupDir = '/work/bioinformatics/shared/dope/export/';
% handles.logfile = 'AnnotationLog.txt'
handles.timeStampStart = char(datetime('now','Format','ddMMMyyyy_hhmm'));
handles.uName = char(java.lang.System.getProperty('user.name'));
handles.compName = char(java.net.InetAddress.getLocalHost.getHostName);
handles.sessionID = [handles.timeStampStart '+' handles.uName '+' handles.compName '_'];
handles.autoSaveCount = 0;
handles.frameUpdatePause = 0.05;
handles.selPtIdx = 1;

% Initialize Label Dictionary
if nargin < 1
%     [FileName,PathName,FilterIndex] = uigetfile(FilterSpec,DialogTitle);
    [filename, pathname] = ...
     uigetfile({'*.mat'},'.MAT CellDataSet File Selector');

    cellDataSet = [pathname filesep filename];
end

if ischar(cellDataSet) && (exist(cellDataSet, 'file') == 2) 

    inM = load(cellDataSet); 
    
    % Capture original file path
    data.info.loadCellDataFile = cellDataSet;
    data.info.sessionID = handles.sessionID;

    % Grab metadata
    if isfield(inM, 'cellDataSet')    
        cellDataSet = inM.cellDataSet;
    elseif isfield(inM, 'allCellsMovieData')    
        cellDataSet = inM.allCellsMovieData;
    end
    cellDataSet = inM.cellDataSet;

    % Load movies... (assumed ot be in same order as meta data)
    if isfield(inM, 'cellMoviesPath') && (exist(inM.cellMoviesPath, 'file') == 2) 

        data.moviesPath = inM.cellMoviesPath;
        choice = questdlg(['Load movies from ' data.moviesPath '?'], ...
                          'Load Movies into Memory?', ...
                          'Yes', ...
                          'No','No');
    else
        choice = 'No';
        
        if isfield(cellDataSet{1},'cellMD') 
            disp('Checking movieData of individual cells...')
            MD = MovieData.load(cellDataSet{1}.cellMD);
            assert(isa(MD, 'MovieData'))
            handles.flyLoad = true;
        end
               
    end 


    if string(choice) == 'Yes'
        hwarn = warndlg(['Please wait...loading ' num2str(length(cellDataSet)) ' movies into memory...']);
        %% TODO detect if MovieList -- if so, OME-TIFF expected... on the fly...
        %% else, give warning that it may take a while to load.
        S = load(data.moviesPath, 'cellMovies');
        if isfield(S, 'cellMovies')    
            data.movies = S.cellMovies;
        elseif isfield(S, 'allCellMovies')    
            data.movies = S.allCellMovies;
        end
        if ishandle(hwarn), close(hwarn), end;
        disp(['Done loading movies.. from ' data.moviesPath]);

    else
        disp('Not loading any images...')
        data.movies = cell(1,length(cellDataSet));
    end

    % pre-defined annotation set (container.Map)
    if isfield(inM, 'annotationSet')
        data.annotationSetIn = inM.annotationSet;
    end

    handles.sessionID = [handles.sessionID '_' data.info.loadCellDataFile];
end

% check if cell array
if iscell(cellDataSet)
    disp('converning to struct array from cell array');
    cellDataSet = cell2mat(cellDataSet);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BoilerPlate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initializeDataStruct_Assaf();
if ~isempty(data.annotationSetIn)
    addAnnotationNoGUI(data.annotationSetIn.keys);
end

% Build out GUI/labels
initMainGUI();

% Append info to UserData
% cxFig = findall(0,'Tag', 'dopeAnnotator');
% cxFig.UserData.handles = handles;
% cxFig.UserData.data = data;
% disp('To update the UserData after events in the GUI, click the "save Notes"');
   
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
    data.meta.anno.set = {}; 
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

    data.extra.time = [cellfun(@(t) t(1), {cellDataSet.ts})]';

    % Initialize backup info....
    handles.backupInfo.cellMoviesPath = data.moviesPath;
    handles.backupInfo.keys = data.meta.key;
    handles.backupInfo.expr = data.meta.class.expr;

    % [MovieData-OMETIFF]
    if handles.flyLoad 
        data.MD = {cellDataSet.cellMD};
        handles.MDcache = cell(1, length(data.meta.mindex));
    end
    
    if string(handles.choiceLoadMD) == 'Yes'
        for i = 1:length(handles.MDcache)
            handles.MDcache{i} = MovieData.load(data.MD{1});
        end
    else
        disp('Not pre-loading any MovieDatas...')
    end 
    
end

%===============================================================================
% Save and Export Functions
%===============================================================================
function writeLog(action, tag, key, expr)

    if exist(handles.logfile, 'file')==2
        fileID = fopen(handles.logfile,'a');
    else
        fileID = fopen('.AnnotationsLogLocalClickFury.txt','a');
    end
    
    if iscell(tag)
        tag = tag{1};
    end
    
    p = mfilename('fullpath');
    [status md5sumout] = system(['md5sum "' p '.m"']);
    timeS = char(datetime('now','Format','ddMMMyyyy hh:mm:ss'));

    formatSpec = '[%s] sessionID[%s] : {%s} \t "%s" \t <%s> \t %s \t %s\n';
   
    fprintf(fileID,formatSpec, timeS, handles.sessionID, action, tag, key, expr, md5sumout);
    fclose(fileID);

    handles.autoSaveCount =  handles.autoSaveCount + 1;
    if (mod(handles.autoSaveCount, 25) == 0)
        exportDataState;
    end
end

%===============================================================================
% Setup main GUI window/figure
%===============================================================================

function initMainGUI()
    
    xsizeF = 1080;
    ysizeF = 720;
      
    % Create main figure
    handles.h1 = figure('Units','pixels', 'Position',[100 100 xsizeF ysizeF],...
                        'Visible',get(0,'defaultfigureVisible'),...
                        'Color',get(0,'defaultfigureColor'),...
                        'CurrentAxesMode','manual',...
                        'IntegerHandle','on',...
                        'MenuBar','none',...
                        'Name','dopeAnnotator',...
                        'NumberTitle','off',...
                        'Tag','dopeAnnotator',...
                        'Resize','off',...
                        'PaperPosition', get(0,'defaultfigurePaperPosition'),...
                        'ScreenPixelsPerInchMode','manual',...
                        'HandleVisibility','callback',...
                        'CloseRequestFcn', @my_closereq);

    function my_closereq(~, ~)
       selection = questdlg({'Close CellXplorer? ','(!) Please verify desired info saved first (!)'},...
          'Close CellXplorer?',...
          'EXIT','RETURN','EXIT'); 
       switch selection
          case 'EXIT'
            warning('add auto save here!')
            delete(gcf)
          case 'RETURN'
          return 
       end
    end
    
    handles.mainP = uipanel('Parent',handles.h1,...
                            'FontUnits','pixels',...
                            'Units','pixels',...
                            'Title','Dope Annotator',...
                            'Tag','mainPanel',...
                            'Position',[5 5 xsizeF-10 ysizeF-15],...
                            'FontSize',12, ...
                            'FontWeight','bold');
    
    [x, y, w, h] = getPosH(handles.mainP);
    
    handles.movie = uipanel('Parent',handles.mainP,...
                              'FontUnits','pixels',...
                              'Units','pixels',...
                              'Title','Cell Movie',...
                              'Tag','uipanel_movie',...
                              'Position',[w-w/2.5-10 h/4 400 400],...
                              'FontSize', 8);
    handles.movie.TitlePosition='righttop';
    
                          
%     [x, y, w, h] = getPosH(handles.mainP);
    hAnnoP = 600;
    wAnnoP = w/2;
    
    handles.annoP = uipanel('Parent',handles.mainP,...
                            'FontUnits','pixels',...
                            'Units','pixels',...
                            'Title','Annotation Panel',...
                            'Tag','uipanel_anno',...
                            'Position',[30 (h-hAnnoP)/2-10 wAnnoP hAnnoP],...
                            'FontSize', 8);
%     handles.annoP.Title='';
%     handles.annoP.BorderType='none';
    
end


%===============================================================================
% Movie display
%===============================================================================

opts = {'Parent', handles.movie, 'Units', 'pixels',...
             'Position',[2 2 handles.movie.Position(3)-2 handles.movie.Position(3)-1],...
             'Color',[1 1 1],'Box' 'off', 'XTick',[],'YTick',[]};

axMovie = axes(opts{:});
axMovie.XColor = 'w';
axMovie.YColor = 'w';
handles.axMovie = axMovie;
set(handles.axMovie, 'XTick', []);
set(handles.axMovie, 'YTick', []);
colormap(handles.axMovie, gray);


function updateMovie()
    MD = handles.MDcache{handles.selPtIdx};
    movieFrame = MD.channels_.loadImage(handles.movies.fidx);

    imagesc(movieFrame,'Parent', handles.axMovie, 'HitTest', 'off');
    set(handles.axMovie, 'XTick', []);
    set(handles.axMovie, 'YTick', []);
    set(handles.axMovie, 'Position', [5 8 size(movieFrame,1)*1.75 size(movieFrame,2)*1.75])    
    colormap(handles.axMovie, gray);
end 


function playMovie(varargin)
    if isempty(handles.MDcache{handles.selPtIdx})
        MD = MovieData.load(data.MD{handles.selPtIdx});
        handles.movies.nf = MD.nFrames_;
        handles.MDcache{handles.selPtIdx} = MD;
    else
        MD = handles.MDcache{handles.selPtIdx};
    end
    handles.movies.nf = MD.nFrames_;
    nf = handles.movies.nf;
    cell_idx = handles.selPtIdx;
    i = 1;
    while (i <= nf) && (handles.selPtIdx == cell_idx) %% ADD SELECTION too...
        handles.movies.fidx = i;
        updateMovie();
        pause(handles.frameUpdatePause);
        i = i+1;
    end
end



function playMovieLoop(varargin)
    tic
    ttoc = 0;
    while ttoc < 5 %% (ADD SELECTION CHECK TOO)
        playMovie;
        ttoc = toc;
    end
    
end

playMovieLoop()
% handles.selPtIdx = 

%===============================================================================
% Annotation Panel Buttons
%===============================================================================

% 
%     %-------------------------------------------------------------------------------
%     % Control/Movie panels of GUI
%     %-------------------------------------------------------------------------------
% 
%     % Start from bottom left
%     % Left to right will be fixed.  
%     % Dynamically expand vertically.
%     % Standard Distance of 5 pixel between panel
%     gapSize = 5;
%     % Data Selection Panel
%     xSizeSelectPanel = 450;
%     % ---> insert dynamic size info here
%     ySizeSelectPanel = 175 + 15*numel(handles.info.labelTypes); 
% 
%     % DR Panel
%     xPosDR = gapSize;
%     yPosDR = ySizeSelectPanel + gapSize*2;
%     xSizeDRPanel= xSizeSelectPanel;
%     ySizeDRPanel = 450;
% 
%     % Movie Panel
%     xSizeLabelPanel = 450;
% 
%     handles.DataSel = uipanel('Parent',handles.mainP,'FontUnits','pixels','Units','pixels',...
%     'Title','Data Selection/View Criterion','Tag','uipanel_select',...
%     'Position',[gapSize gapSize xSizeSelectPanel ySizeSelectPanel],...
%     'FontSize',13);
% 
%     handles.h2_DR = uipanel('Parent',handles.mainP, 'FontUnits','pixels', 'Units','pixels',...
%     'Title','2D Visualization - Dimension Reduction',...
%     'Tag','uipanel_axes',...
%     'Position',[handles.DataSel.Position(1), handles.DataSel.Position(4)+handles.DataSel.Position(2)+gapSize, xSizeDRPanel, ySizeDRPanel],...
%     'FontSize',13);
% 
%     handles.LabelA = uipanel('Parent',handles.mainP,'FontUnits','pixels','Units','pixels',...
%     'Title','Cell Labeling',...
%     'Tag','uipanel_annotate',...
%     'Position',[handles.h2_DR.Position(3)+handles.h2_DR.Position(1)+gapSize,...
%                 handles.DataSel.Position(2),...
%                 xSizeLabelPanel 475],...
%     'FontSize',13);
% 
% 
%     widthCellInfo = 330;
%     handles.widthCellInfo = widthCellInfo;
% 
%     xposLabels = widthCellInfo/2;
%     LabelH = 13;
%     gapL = 3;
%     widthString = 55;
%     handles.cellInfo = uipanel(...
%     'Parent',handles.mainP,...
%     'FontUnits','pixels',...
%     'Units','pixels',...
%     'Title','Cell Info',...
%     'Tag','CellInfopanel',...
%     'Position',[handles.h2_DR.Position(3)+handles.h2_DR.Position(1)+gapSize,...
%                 handles.LabelA.Position(2)+handles.LabelA.Position(4)+gapSize,...
%                 handles.LabelA.Position(3),...
%                 handles.mainP.Position(4)-handles.LabelA.Position(4)-30],...
%     'FontSize',13);
% 
% 

% 
% 
 
%     
%     
%     % reconstitute original cell array with annotations
% %     function re
%     
%     
% %-------------------------------------------
% % 
% % 
%     handles.ActionNotice = uicontrol(...
%     'Parent',handles.mainP,...
%     'FontUnits','pixels',...
%     'Units','pixels',...
%     'Style','text',...
%     'String','Please wait updating Plots...',...
%     'Position',[handles.mainP.Position(3)-250 100 250 35],...
%     'Tag','SaveNotes',...
%     'FontSize',15,...
%     'BackgroundColor',[.75 .5 1]);
%     set(handles.ActionNotice, 'Visible', 'off');

%===============================================================================
% Annotation Management functions
%===============================================================================   

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
            writeLog('create-tag', newStr, 'none', 'none');
        end
    end
end

% % %===============================================================================
% % % Helper functions
% % %===============================================================================   

function [x, y, w, h]=getPosH(Hin)
   x = Hin.Position(1);
   y = Hin.Position(2);
   w = Hin.Position(3); 
   h = Hin.Position(4); 
end

end