function [cxFig handles] = dopeAnnotator(cellDataSet, varargin)
%DOPEANNOTATOR rapid GUI annotation labeling tool for 2D cell movies
%
% Inputs:
% 		   cellDataSet:	struct or .mat file containing struct with cellMovie metadata 
%
% Andrew R. Jamieson, Dec. 2016 updated JUNE 2017
%
ff = findall(0,'Tag','dopeAnnotator');delete(ff)
data.annotationSetIn = {};

% actual cell sequence (grows as user annotatates) (for this session)
data.cellsAnnotatedList = [];

if nargin >= 1
    ip = inputParser;
    ip.KeepUnmatched = true;
    ip.CaseSensitive = false;
    ip.addRequired('cellDataSet', @(x) isstruct(x) || iscell(x) || exist(x, 'file')==2 || @(x) isa(x, 'MovieList'));
    ip.addOptional('cellAnnotations', {}, @(x) isstruct(x) || exist(x, 'file')==2);
    
    ip.addParameter('annotationSet', {}, @(x) isa(x,'containers.Map'));
    
    ip.parse(cellDataSet, varargin{:});
    data.cellAnnotationsPrev = ip.Results.cellAnnotations;
    data.annotationSetIn = ip.Results.annotationSet;
   
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
handles.backupDir = '/work/bioinformatics/shared/dope/export/';
handles.logfile = '/work/bioinformatics/shared/dope/export/.AnnotationsLog.txt';
handles.timeStampStart = char(datetime('now','Format','ddMMMyyyy_hhmm'));
handles.uName = char(java.lang.System.getProperty('user.name'));
handles.compName = char(java.net.InetAddress.getLocalHost.getHostName);
handles.sessionID = [handles.timeStampStart '-' handles.uName '-' handles.compName '_'];
handles.sessionID2 = [handles.timeStampStart '-' handles.uName];
handles.Masterlogfile = ['/work/bioinformatics/shared/dope/export/MasterAnnotationsLogClickFury_' handles.sessionID2 '.txt'];
handles.matlabSaveFile = ['/work/bioinformatics/shared/dope/export/backup_dopeAnnotator_output_' handles.sessionID2 '.mat'];
handles.autoSaveCount = 0;
handles.frameUpdatePause = 0.075;
handles.movieLoopLimit = Inf;
handles.selPtIdx = 1;
handles.stageDriftCorrection = true;
handles.maxPerRow = 4;
handles.numRows = 2;
handles.buttonSizeH = 150;
handles.buttonSizeW = 150;
handles.repeatsAllowed = true;
handles.pauseProgress = false;


% Initialize Label Dictionary
if nargin < 1
    [filename, pathname] = ...
     uigetfile({'*.mat'},'.MAT CellDataSet File Selector');
    cellDataSet = [pathname filesep filename];
end

if ischar(cellDataSet) && (exist(cellDataSet, 'file') == 2) 

    inM = load(cellDataSet); 
    % Capture original file path
    data.info.loadCellDataFile = cellDataSet;
    data.info.sessionID = handles.sessionID;

    % Grab main metadata
    cellDataSet = inM.cellDataSet;
    
    if isfield(cellDataSet{1},'cellMD') 
        disp('Checking movieData of individual cells...')
        MD = MovieData.load(cellDataSet{1}.cellMD);
        assert(isa(MD, 'MovieData'))
        handles.flyLoad = true;
    end               

    % pre-defined annotation set (container.Map)
    if isfield(inM, 'annotationSet')
        data.annotationSetIn = inM.annotationSet;
    end

    handles.sessionID = [handles.sessionID '_' data.info.loadCellDataFile];
else
    disp('Using variable passed cellDataSet to dopeAnnotator');
end

% check if cell array
if iscell(cellDataSet)
    disp('converning to struct array from cell array');
    cellDataSet = cell2mat(cellDataSet);
end

% check if previous annotations
if isempty(data.cellAnnotationsPrev)
    choice = questdlg(['Load previous annotation set?'],...
                       'Load previous annotation set?',...
                       'Yes','No','No');
    if string(choice) == 'Yes'
        [filename, pathname] = ...
         uigetfile({'*.mat'},'.MAT cellAnnotations File Selector');
        data.cellAnnotationsPrev = [pathname filesep filename];
    end
end

data.remainingKeys = [];
if ~isempty(data.cellAnnotationsPrev)
    if ischar(data.cellAnnotationsPrev) && (exist(data.cellAnnotationsPrev, 'file') == 2) 
        inA = load(data.cellAnnotationsPrev);
        if isfield(inA, 'cellAnnotations')    
            data.cellAnnotationsPrev = inA.cellAnnotations;
            if ~isfield(data.cellAnnotationsPrev, 'sessionIDList')
                data.cellAnnotationsPrev.sessionIDList = {'-undefined-'};
                warning('No previous sessionIDList found...adding empty {-undefined-}');
            end
        else
            error('Expected structure not found: cellAnnotations');
        end    
    end
        
    disp('Truncating cellDataSet based on previous...annotations');
    doneAnnotations = data.cellAnnotationsPrev;
    tempDoneKeys = {doneAnnotations.cellData.key};
    % check for empty annotations (i.e, "timeouts") to exclude
    indx_tO = [doneAnnotations.cellData.timeOutFlag] == 1;
    doneKeys = tempDoneKeys(~indx_tO);
    disp(['# of previously COMPLETED keys : ' num2str(length(doneKeys))]);
    
    % find remaining keys
    allCellKeys = {cellDataSet.key}; % from total dataset
    disp(['# of ORIGINAL data set keys: ' num2str(length(unique(allCellKeys)))]);
    remainingKeys = setxor(allCellKeys, doneKeys);
    disp(['# of REMAINING keys : ' num2str(length(remainingKeys))]);
    data.remainingKeys = remainingKeys;
    subSetData = cell(1, length(remainingKeys));

    % create SubSet
    j = 1;
    for i = 1:length(cellDataSet)
        if ismember({cellDataSet(i).key}, remainingKeys)
            subSetData{j} = cellDataSet(i);
            j = j + 1;
        else
            data.cellsAnnotatedList = [data.cellsAnnotatedList i];
        end
    end

    % cellDataSet = subSetData;

    handles.info.flagPrevAnnotations = true;
    handles.info.prevAnnotationSessionIDs = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ask User for annotation configs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prompt for # of cells to annotation (out of total)
% prompt for % of repeats to annotate.
function promptAnnoConfig(varargin)
    defaultans = {num2str(length(cellDataSet)),'10'};
    goodAns = false;
    while ~goodAns
        ans = inputdlg({['How many annoations to conduct? (out of ' num2str(length(cellDataSet)) ...
            ' Total Cells)'], 'probability of repeats? enter [0-100]'},...
                      'Configure Annotatotion Sequence', [1 50; 1 50], defaultans);
        if (str2num(ans{2}) <= 100) && (str2num(ans{2}) >= 0) && (str2num(ans{1}) > 0) &&  (str2num(ans{1}) <= length(cellDataSet))
            goodAns = true;
        end
        handles.prctRepeats = str2num(ans{2}); 
        handles.numCellsToAnnotate = str2num(ans{1});    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BoilerPlate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
handles.repeatFlag = 0;
handles.junkFlag = 0;
handles.timeOutFlag = 0;
handles.firstSelectionMade = 0;
handles.NextCell = 0;

initializeDataStruct_Assaf();
promptAnnoConfig()
configAnnoSequence();
if ~isempty(data.annotationSetIn)
    addAnnotationNoGUI(data.annotationSetIn.keys);
end
% Build out GUI/labels
initMainGUI();
movieInit();
presentCells;

% SAVE
% CLOSE ?
   
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
    % These will model behaviior of original CellX unique for each cell
    % 
    % Note this is different than the simple sequence of annotations 
    % where repeats can occur and different mappings.
    data.meta.anno.tagMap = containers.Map('KeyType','char','ValueType', 'any'); % tags to cells (by local master index)
    data.meta.anno.tagMapKey = containers.Map('KeyType','char','ValueType', 'any'); % tags to cells(by unique Key)
    data.meta.anno.RevTagMap = containers.Map('KeyType','char','ValueType', 'any'); % cell to tags
    % Set specific to each cell in array
    data.meta.anno.byCell = cell(1,length(data.meta.key)); % will be generated at export or save 

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

function configAnnoSequence(varargin)

    handles.prctRepeats;
    handles.numCellsToAnnotate;
    data.seedCellSeq = randsample(length(data.meta.mindex), length(data.meta.mindex), false);
         
    % This will store the annotations with a struct containing:    
    data.cellAnnotationsData.cellData = [];
    data.cellAnnotationsDataDetail = [];
    data.cellAnnotationsData.sessionIDList = {handles.sessionID};
    data.cellAnnotationsData.Start_cellsAnnotatedList = data.cellsAnnotatedList;

    if ~isempty(data.cellsAnnotatedList)
        data.cellAnnotationsData.cellData = [data.cellAnnotationsPrev.cellData];
        data.cellAnnotationsData.sessionIDList = {data.cellAnnotationsPrev.sessionIDList{:} handles.sessionID};
    end

end 

%===============================================================================
% Save and Export Functions
%===============================================================================
function saveMatFile(varargin)
    cellAnnotations = data.cellAnnotationsData;
    cellAnnotations.End_cellsAnnotatedList = data.cellsAnnotatedList;
    file1_name = ['dopeAnnotator_output_' handles.sessionID2 '.mat'];
    save(file1_name, 'cellAnnotations');
    save(handles.matlabSaveFile, 'cellAnnotations');
    if (mod(handles.autoSaveCount, 99) == 0)
        handles.ActionNotice.String = '{SAVING}';
        handles.ActionNotice.BackgroundColor = [0 0 1];
        set(handles.ActionNotice, 'Visible', 'on');
        cellAnnotations2 = data.cellAnnotationsDataDetail;
        save(['Detail_dopeAnnotator_output_' handles.sessionID2 '.mat'], 'cellAnnotations2');
        save(handles.matlabSaveFile, 'cellAnnotations');
        set(handles.ActionNotice, 'Visible', 'off');
    end      
    disp(['done ...saving...' handles.matlabSaveFile]);
end

function writeMasterLog(action, tag, key, expr)
    % if exist(handles.logfile, 'file')==2
    %     fileID = fopen(handles.Masterlogfile,'a');
    fileID2 = fopen(['MasterAnnotationsLogClickFury_' handles.timeStampStart '.txt'],'a');
    fileID3 = fopen(handles.logfile,'a');
    % end
    loopTags = false;
    skipCellstr = num2str(0);
    if iscell(tag) && length(tag) == 1
        tag = tag{1};
    elseif iscell(tag) && length(tag) > 1
        loopTags = true;
        tags = tag;
    elseif isempty(tag) && iscell(tag)
        tag = 'NONE';
        if handles.NextCell == 1
            skipCellstr = num2str(handles.NextCell);
        else
            skipCellstr = num2str(0);
        end
    end
    
    p = mfilename('fullpath');
    [status md5sumout] = system(['md5sum "' p '.m"']);
    timeS = char(datetime('now','Format','ddMMMyyyy hh:mm:ss'));
   
    formatSpec = 'dope[%s]\t[sessionID][%s]\t{%s} \t"%s"\tskipCell:%s\ttimeoutFlag:%s\tjunkFlag:%s\trepeatFlag:%s\t<%s>\t%s\t%s\n';
   
    if loopTags
        for iT = 1:length(tags)
    %             fprintf(fileID, formatSpec, timeS, handles.sessionID, action, tags{iT},skipCellstr, num2str(handles.timeOutFlag),num2str(handles.junkFlag),num2str(handles.repeatFlag), key, expr, md5sumout);
            fprintf(fileID2, formatSpec, timeS, handles.sessionID, action, tags{iT},skipCellstr, num2str(handles.timeOutFlag),num2str(handles.junkFlag),num2str(handles.repeatFlag), key, expr, md5sumout);
            fprintf(fileID3, formatSpec, timeS, handles.sessionID, action, tags{iT},skipCellstr, num2str(handles.timeOutFlag),num2str(handles.junkFlag),num2str(handles.repeatFlag), key, expr, md5sumout);
        end
    else
    %         fprintf(fileID, formatSpec, timeS, handles.sessionID, action, tag, skipCellstr, num2str(handles.timeOutFlag),num2str(handles.junkFlag),  num2str(handles.repeatFlag), key, expr, md5sumout);
        fprintf(fileID2, formatSpec, timeS, handles.sessionID, action, tag, skipCellstr, num2str(handles.timeOutFlag),num2str(handles.junkFlag),  num2str(handles.repeatFlag), key, expr, md5sumout);
        fprintf(fileID3, formatSpec, timeS, handles.sessionID, action, tag, skipCellstr, num2str(handles.timeOutFlag),num2str(handles.junkFlag),  num2str(handles.repeatFlag), key, expr, md5sumout);
    end
    
    %     fclose(fileID);
    fclose(fileID2);
    fclose(fileID3);
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
      saveMatFile; 
        switch selection
          case 'EXIT'
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

    handles.cellMDtext = uicontrol(...
    'Parent',handles.mainP,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'HorizontalAlignment','right',...
    'String', ['Cell MovieData :' data.MD{handles.selPtIdx}], ...
    'Style','text',...
    'Position',[handles.mainP.Position(3)-1000 handles.mainP.Position(4)-35, 950, 15],...
    'Tag','cellMDtext',...
    'FontSize',8);
    set(handles.cellMDtext, 'Visible', 'off');

    handles.masterMDtext = uicontrol(...
    'Parent', handles.mainP,...
    'FontUnits','pixels',...
    'Units','pixels',...
    'HorizontalAlignment','right',...
    'String', ['Master MovieData :' data.meta.class.MD{handles.selPtIdx}], ...
    'Style','text',...
    'Position',[handles.mainP.Position(3)-1000 handles.mainP.Position(4)-45, 950, 15],...
    'Tag','masterMDtext',...
    'FontSize',8);
    set(handles.masterMDtext, 'Visible', 'off');

    
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
    handles.annoP.Title='';
    handles.annoP.BorderType='none';

    %===============================================================================
    % Annotation Panel Buttons
    %===============================================================================
    % -------------------------------------------------------------------------------    

    bH = wAnnoP/4-1;%handles.buttonSizeH;
    bW = wAnnoP/4-1;% handles.buttonSizeW;
    
    function NextButton(varargin)
        handles.NextCell = 1;
    end

    buttonOptsNext = {'Parent', handles.mainP,...
                  'FontUnits','pixels',...
                  'Units','pixels',...
                  'HorizontalAlignment','center',...
                  'FontSize',22,...
                  'Style','pushbutton',...
                  'BackgroundColor', [.91 .93 .93],...
                  'Callback',@NextButton};              
              
    handles.nextButton = uicontrol(buttonOptsNext{:},...
                                  'String','NEXT CELL',...
                                  'Position',[w/2+100 2 bW bH],...
                                  'Tag','nextButtonTag');
    
    function buttonSelected(src, ~)        
        if src.Value == 1
            src.BackgroundColor = [1 0 0];
%             src.ForegroundColor = [1 1 0];
            src.ForegroundColor = [0 1 1];
%             src.ForegroundColor = [0 0 0];
        elseif src.Value == 0
            src.ForegroundColor = [0 0 0];
            src.BackgroundColor = [.87 .84 .84];
        end
        
%         cellKey = data.meta.key{handles.selPtIdx};
%         tagDataPointNoGUI(src.String, handles.selPtIdx, cellKey, src.Value);
        
        handles.firstSelectionMade = 1;
        
        if strcmp(src.String, handles.bFocus.String) || ...
                strcmp(src.String, handles.bJunk.String) || ...
                strcmp(src.String,  handles.bUD.String) || ...
                strcmp(src.String, handles.bMulti.String)
        
            handles.junkFlag = 1;
        
        end
        
        updateStatus;
    end
    
    % -------------------------------------------------------------------------------
    % Create Standard Buttons (always present)
    % -------------------------------------------------------------------------------

    yPos = 25;
    buttonOpts = {'Parent', handles.annoP,...
                  'FontUnits','pixels',...
                  'Units','pixels',...
                  'HorizontalAlignment','center',...
                  'FontSize',22,'Style','togglebutton',...
                  'BackgroundColor', [.87 .84 .84],...
                  'Value', 0,...
                  'Callback',@buttonSelected};

              
              
    handles.bJunk = uicontrol(buttonOpts{:},...
                                  'String','JUNK',...
                                  'Position',[1+bW*0 yPos bW bH],...
                                  'Tag','junkButtonDefault');

    handles.bUD = uicontrol(buttonOpts{:},...
                                'String','UNDEFINED',...
                                'Position',[1+bW*1+1 yPos bW bH],...
                                'Tag','undefinedButtonDefault');

    handles.bFocus = uicontrol(buttonOpts{:},...
                                'Position',[1+bW*2+1 yPos bW bH],...
                                'Tag','focusButtonDefault');
    set(handles.bFocus,'String','<html>OUT<br>OF<br>FOCUS');

    handles.bMulti = uicontrol(buttonOpts{:},...
                              'Position',[1+bW*3+1 yPos bW bH],...
                              'Tag','multiCellButtonDefault');    
    set(handles.bMulti,'String','<html>MULTIPLE<br>CELLS');

    % -------------------------------------------------------------------------------
    % Create dynamic Annotation labels
    % -------------------------------------------------------------------------------   
    buttonOpts2 = {'Parent', handles.annoP,...
                  'FontUnits','pixels',...
                  'Units','pixels',...
                  'HorizontalAlignment','center',...
                  'FontSize',20,'Style','togglebutton',...
                  'BackgroundColor', [.92 .92 .92], ...
                  'Value', 0,...
                  'Callback',@buttonSelected};    
   
    numA = numel(data.meta.anno.set);
    handles.annoB = gobjects([numA, 1]);
    
    yPos = yPos + bH + 15;
    xPos = 0;
    iA = 1;
    while iA <= numA
        
        handles.annoB(iA) = uicontrol(buttonOpts2{:},...
                                      'String', data.meta.anno.set{iA},...
                                      'Position',[xPos+1 yPos bW bH],...
                                      'Tag',[data.meta.anno.set{iA} 'DynamicTag']);
                
      xPos = xPos + bW;
      if mod(iA,4) == 0
            yPos = yPos + bH;
            xPos = 0;
      end
      iA = iA + 1;
    end
 
%     %-------------------------------------------------------------------------------
%     % Control/Movie panels of GUI
%     %-------------------------------------------------------------------------------
% 
    [x, y, w, h] = getPosH(handles.mainP);
    handles.ActionNotice = uicontrol(...
                                    'Parent',handles.mainP,...
                                    'FontUnits','pixels',...
                                    'Units','pixels',...
                                    'Style','text',...
                                    'String','All done!...',...
                                    'Position',[w/2+250 60 250 55],...
                                    'Tag','NoticeWarnText',...
                                    'FontSize',50,...
                                    'BackgroundColor',[.75 .5 1]);
    set(handles.ActionNotice, 'Visible', 'off');
   
    handles.annoBDefaults = findall(0,'-regexp','Tag', 'ButtonDefault');
    handles.annoBDynam    = findall(0,'-regexp','Tag', 'Dynamic');
end

%===============================================================================
% Cell Array Management
%===============================================================================

    function presentCells(varargin)
        numAnnotated = 0;
        iOrigSeq = 1;
        while numAnnotated < handles.numCellsToAnnotate
            resetAnnotations();
            if handles.repeatsAllowed && ~isempty(data.cellsAnnotatedList) && rand <= handles.prctRepeats/100
                newCell = randsample(unique(data.cellsAnnotatedList), 1, false);
                handles.repeatFlag = 1;
            else
                newCell = data.seedCellSeq(iOrigSeq);
                iOrigSeq = iOrigSeq + 1;
            end

            handles.selPtIdx = newCell;
            handles.cellMDtext.String = ['Cell MovieData :' data.MD{handles.selPtIdx}];
            handles.masterMDtext.String = ['Master MovieData :' data.meta.class.MD{handles.selPtIdx}];
            
            playMovieLoop;
            if (handles.firstSelectionMade == 0) && (handles.NextCell == 0)
                handles.timeOutFlag = 1;
            end
            data.cellsAnnotatedList = [data.cellsAnnotatedList newCell];
            numAnnotated = numAnnotated + 1;
            collectAnnotations(newCell);

            handles.autoSaveCount =  handles.autoSaveCount + 1;
            if (mod(handles.autoSaveCount, 10) == 0)
                handles.ActionNotice.String = '{SAVING}';
                handles.ActionNotice.BackgroundColor = [0 0 1];
                set(handles.ActionNotice, 'Visible', 'on');
                saveMatFile;
                set(handles.ActionNotice, 'Visible', 'off');
            end            
            
            % write log
            % save file 
        end
        handles.ActionNotice.String = 'All DONE!';
        handles.ActionNotice.BackgroundColor = [0 1 0];
        set(handles.ActionNotice, 'Visible', 'on');
        
        % SAVE 
        % CLOSE!
    end

    function addCellAnnotation(cell_index, annotations)
        cellKey = data.meta.key{cell_index};
        cellexpr = data.meta.expStr{cell_index};
%         assert(strcmp(cellKey,data.seedCellSeq_byKey(cell_index)))
        newCell = struct();
        
        newCell.key = cellKey;
        newCell.repeatFlag = handles.repeatFlag;
        newCell.junkFlag = handles.junkFlag;
        newCell.timeOutFlag = handles.timeOutFlag;
        newCell.NextCellFlad = handles.NextCell;
        newCell.annotations = annotations;
        newCell.cellexpr = cellexpr;

        if  length(data.meta.anno.RevTagMap.keys) > 0
            newCell.annotationSet = data.meta.anno.RevTagMap(cellKey);
        else
            newCell.annotationSet = {};
        end
        newCell2 = newCell;
        newCell2.cellMD = data.MD{cell_index};
        newCell2.sessionID = handles.sessionID;

        data.cellAnnotationsDataDetail = [data.cellAnnotationsDataDetail newCell2];
        data.cellAnnotationsData.cellData = [data.cellAnnotationsData.cellData newCell];
        writeMasterLog('add', annotations, cellKey, cellexpr);
    end

    function updateCellMD(cell_index, annotations)
        cellKey = data.meta.key{handles.selPtIdx};
        MD = MovieData.loadMatFile(data.MD{cell_index});
        extProcName = ['DopeAnnotations_' handles.timeStampStart];

        extProcParams.cellKey = cellKey;
        extProcParams.annotations = annotations;
        extProcParams.sessionID = handles.sessionID;
        extProcParams.timeStampStart = handles.timeStampStart;
        
        if isempty(MD.getProcessIndex(extProcName))
            extProc = ExternalProcess(MD, extProcName);
            extProc.setParameters(extProcParams);
            MD.addProcess(extProc);
            MD.save();
        else
            extProcindx = MD.getProcessIndex(extProcName);
            MD.processes_{extProcindx}.setParameters(extProcParams);
            MD.save();
        end
    end

%===============================================================================
% Annotation Management functions
%===============================================================================   

    function resetAnnotations(varargin)
                      
        for i = 1:length(handles.annoBDefaults)
            handles.annoBDefaults(i).Value = 0;
            handles.annoBDefaults(i).ForegroundColor = [0 0 0];
            handles.annoBDefaults(i).BackgroundColor = [.87 .84 .84];    
        end
        
        for i = 1:length(handles.annoBDynam)
            handles.annoBDynam(i).Value = 0;
            handles.annoBDynam(i).ForegroundColor = [0 0 0];
            handles.annoBDynam(i).BackgroundColor = [.92 .92 .92];
        end
        handles.repeatFlag = 0;
        handles.junkFlag = 0;
        handles.timeOutFlag = 0;
        handles.firstSelectionMade = 0;
        handles.NextCell = 0;
    end

    function collectAnnotations(cell_index)
        annotations = {};
        for i = 1:length(handles.annoBDefaults)
            tbut = handles.annoBDefaults(i);
            if tbut.Value
                annotations = [annotations {tbut.String}];
                handles.junkFlag = 1;
            end
        end
        
        for i = 1:length(handles.annoBDynam)
            tbut = handles.annoBDynam(i);
            if tbut.Value
                annotations = [annotations {tbut.String}];
            end
        end
        addCellAnnotation(cell_index, annotations)
        updateCellMD(cell_index, annotations)
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
            writeMasterLog('create-tag', newStr, 'none', 'none');
        end
    end
end

%===============================================================================
% Movie display
%===============================================================================

function movieInit(varargin)
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


    [x, y, w, h] = getPosH(handles.mainP);
    

    handles.toggleSDC = uibuttongroup('Parent',handles.mainP,...
    'FontUnits','pixels','Units','pixels',...
    'Title','Stage Drift Correction',...
    'Tag','toggleSDC',...
    'Position',[w-w/2.5-10 h/5 150 35], ...
    'FontSize', 12,...
    'SelectionChangedFcn', @sdcCallback);    

    handles.SDC1 = uicontrol(handles.toggleSDC,'Style','radiobutton',...
                      'FontUnits','pixels','Units','pixels',...
                      'String','on',...
                      'Position',[15 1 35 20],...
                      'HandleVisibility','off',...
                      'FontSize', 12);

    handles.SDC0 = uicontrol(handles.toggleSDC,'Style','radiobutton',...
                      'FontUnits','pixels','Units','pixels',...
                      'String','off',...
                      'Position',[50 1 35 20],...
                      'HandleVisibility','off','FontSize', 12);

    % Default to the OR switch    
    handles.toggleSDC.SelectedObject = handles.SDC1;

    function sdcCallback(varargin)
        if strcmp(handles.toggleSDC.SelectedObject.String,'on')
            handles.stageDriftCorrection = true; 
        else
            handles.stageDriftCorrection = false;
        end

    end
    handles.toggleSDC.Visible = 'off';

    rates = string([1 .5 .25 .1 .075 .05 .025 .01 .005]);
    handles.controlFrameRate = uicontrol(...
                                        'Parent',handles.mainP,...
                                        'FontUnits','points',...
                                        'FontSize',8,...
                                        'String', rates, ...
                                        'Style','popupmenu',...
                                        'Value', find(ismember(rates, string(handles.frameUpdatePause))),...
                                        'Units','pixels',...
                                        'Tag','pointSize',...
                                        'Position',[(w-w/2)+370 h/5+5 50 25],...
                                        'Callback',@updateFrameRate);

    handles.frText = uicontrol('Parent',handles.mainP,...
                       'FontUnits','points',...
                       'FontSize',8,...
                       'String', 'Frame Rate Pause', ...
                       'Style','text',...
                       'Position',[(w-w/2)+260 h/5 110 25]);


    function updateFrameRate(source, ~)
       val = source.Value;
       maps = source.String;
       handles.frameUpdatePause = str2double(maps{val});
       disp(['Updating frame pause to : ', (maps{val})]);
    end

    handles.controlFrameRate.Visible = 'off';
    handles.frText.Visible = 'off';

    
    loopTime = string([1 2 3 4 5 6 7 10 Inf]);
    handles.controlLoopTime = uicontrol(...
                                        'Parent',handles.mainP,...
                                        'FontUnits','points',...
                                        'FontSize',8,...
                                        'String', loopTime, ...
                                        'Style','popupmenu',...
                                        'Value', find(ismember(loopTime, string(handles.movieLoopLimit))),...
                                        'Units','pixels',...
                                        'Tag','pointSize',...
                                        'Position',[(w-w/2)+370 h-103 50 25],...
                                        'Callback',@updateLoopLimit);

    handles.loopTimeText = uicontrol('Parent',handles.mainP,...
                       'FontUnits','points',...
                       'FontSize',8,...
                       'String', 'Loop time limit', ...
                       'Style','text',...
                       'Position',[(w-w/2)+275 h-100 95 20]);


    function updateLoopLimit(source, ~)
       val = source.Value;
       maps = source.String;
       handles.movieLoopLimit = str2double(maps{val});
       disp(['Updating loop limit to : ', (maps{val})]);
    end
    
    
    handles.controlLoopTime.Visible = 'off';
    handles.loopTimeText.Visible = 'off';
    
    set(handles.h1, 'KeyPressFcn', {@pb_fig, handles.h1});
    set(handles.annoP.Children, 'KeyPressFcn', {@pb_fig, handles.h1});

    function pb_fig(varargin)
        switch varargin{1,2}.Character
            case 's'
                if strcmp(handles.toggleSDC.Visible,'off')
                    handles.toggleSDC.Visible = 'on';
                else
                    handles.toggleSDC.Visible = 'off';
                end
            case  'p'
                if handles.pauseProgress
                    handles.pauseProgress = false;
                    set(handles.ActionNotice, 'Visible', 'off');
                else
                    handles.pauseProgress = true;
                    set(handles.ActionNotice, 'String', '[!] PAUSED [!]');
                    handles.ActionNotice.BackgroundColor = [1 0 0];
                    set(handles.ActionNotice, 'Visible', 'on');
                end
            case 't'
                if strcmp(handles.controlLoopTime.Visible,'off')
                    handles.controlLoopTime.Visible = 'on';
                    handles.loopTimeText.Visible = 'on';
                else
                    handles.controlLoopTime.Visible = 'off';
                    handles.loopTimeText.Visible = 'off';
                end
            case 'f'
                if strcmp(handles.controlFrameRate.Visible,'off')
                    handles.controlFrameRate.Visible = 'on';
                    handles.frText.Visible = 'on';
                else
                    handles.controlFrameRate.Visible = 'off';
                    handles.frText.Visible = 'off';
                end
            case 'n'
                handles.NextCell = 1;
            case 'l'
                if strcmp(handles.cellMDtext.Visible, 'off')
                    handles.cellMDtext.Visible = 'on';
                    handles.masterMDtext.Visible = 'on';
                else
                    handles.cellMDtext.Visible = 'off';
                    handles.masterMDtext.Visible = 'off';
                end                
            case 'e'
                disp('Saving...');
                assignin('base', 'handlesCX', handles);
                assignin('base', 'dataCX', data);
                saveMatFile;
            otherwise
                disp('button pressed but not utilized');
        end
    end
    
end

function updateMovie()

    MD = handles.MDcache{handles.selPtIdx};
    if handles.stageDriftCorrection 
        SDCindx = MD.getProcessIndex('EfficientSubpixelRegistrationProcess');
        if ~isempty(SDCindx)
            movieFrame = MD.processes_{SDCindx}.loadOutImage(1, handles.movies.fidx);
        else
            movieFrame = MD.channels_.loadImage(handles.movies.fidx);
            handles.toggleSDC.SelectedObject = handles.SDC0;
            handles.stageDriftCorrection = false;
            warning('no Stage Drift Correction Process found: disabling');
        end
    else
        movieFrame = MD.channels_.loadImage(handles.movies.fidx);
    end
    % if mod(handles.selPtIdx, 2)
    %     assert(false)
    % end
    
    imagesc(movieFrame,'Parent', handles.axMovie, 'HitTest', 'off');
    set(handles.axMovie, 'XTick', []);
    set(handles.axMovie, 'YTick', []);
    set(handles.axMovie, 'Position', [5 8 size(movieFrame,1)*1.75 size(movieFrame,2)*1.75])    
    colormap(handles.axMovie, gray);
end 

function playMovie(varargin)
    if isempty(handles.MDcache{handles.selPtIdx})
        MD = MovieData.loadMatFile(data.MD{handles.selPtIdx});
        handles.movies.nf = MD.nFrames_;
        handles.MDcache{handles.selPtIdx} = MD;
    else
        MD = handles.MDcache{handles.selPtIdx};
    end
    handles.movies.nf = MD.nFrames_;
    nf = handles.movies.nf;
    cell_idx = handles.selPtIdx;
    i = 1;
    while (handles.junkFlag == 0) && (i <= nf) && (handles.selPtIdx == cell_idx) && (handles.ttoc < handles.movieLoopLimit) && ...
        (handles.NextCell == 0)
        handles.movies.fidx = i;
        try 
            updateMovie();
        catch
            warning('MD failure!');
            MD
            cell_idx
            handles.movies.fidx
            warning('Movie fail to load properly!');
            handles.pauseProgress = true;
            set(handles.ActionNotice, 'String', '[!] PAUSED [!]');
            handles.ActionNotice.BackgroundColor = [1 0 0];
            set(handles.ActionNotice, 'Visible', 'on');
            uiwait(msgbox({['[ERROR with this movie]: ' MD.movieDataPath_] ' ' 'record MD info & tell ARJ' ' ' 'click NEXT CELL then close this message when ready to proceed' 'close  mmesage'}));
        end

        pause(handles.frameUpdatePause);
        if handles.pauseProgress
            pause(.25);
            uiwait(msgbox({'program paused - ' 'click on red [!] and then keypress p to resume'  '(and then close this box)'}));
        else
            i = i + 1; 
        end
    end
end

function playMovieLoop(varargin)
    tic;
    handles.ttoc = 0;
    cell_idx = handles.selPtIdx;
    while (handles.ttoc < handles.movieLoopLimit) && (handles.selPtIdx == cell_idx) && (handles.junkFlag == 0) ...
            && (handles.NextCell == 0)
        playMovie;
        handles.ttoc = toc;
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

    function updateStatus(varargin)
        % Append info to UserData
        cxFig = findall(0,'Tag', 'dopeAnnotator');
        cxFig.UserData.handles = handles;
        cxFig.UserData.data = data;
    end
end