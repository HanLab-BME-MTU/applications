function label_loadmovieCB(varargin);
%loads movie to show in labelgui

%find position of previous panels
contFigH = findall(0,'Tag','labelgui'); %controlFigure = labelgui
positions = GetUserData(contFigH,'positions'); %get position struct
%look for existing labelPanel first
labelPanelH = findall(0,'Tag','LabelPanel');
if ~isempty(labelPanelH)
    labelPanelPos = get(labelPanelH(end),'Position');
elseif isfield(positions,'labelPanelPos')
    labelPanelPos = positions.labelPanelPos;
else
    labelPanelPos = [];
end



if nargin==0|isempty(varargin{1}) %no (load-all&exist projectData)
    %check if default biodata-dir exists and cd if exist
    mainDir = cdBiodata(2);
    
    %load movie file
    [fname,fpath] = uigetfile({'*.r3d;*.r3c;moviedat*;*.fim',  'All moviefiles'},'select image file');
    if(fname(1)==0)
        image = [];
        fname = [];
        return;
    end;
    cd(fpath);
    filename = [fpath fname];
    
    if isempty(findstr(fname,'.r3d'))
        [mv,stat] = readmat(filename); %works for *.fim and *.r3c and moviedat
    else
        [mv] = r3dread(filename);
    end;
else %filtered movie is being passed down
    mv = varargin{1};
    fname = varargin{2};
    fpath = varargin{3};
end

if ~isempty(mv)
    
    % store movie
    
    % display movie
    % intialize params
    iptsetpref('ImshowBorder', 'tight');
    
    %make movie span the range [0,1]
    %problem: there can be one (or several) very bright frames that should
    %be excluded. Hence take maxima over each frame and exclude everything
    %that is more than 4 std over the mean intensity
    
    %do minimum subtraction
    mv = mv-min(mv(:));
    
    %calculate maximum
    movieSize = size(mv);
    reshMovie = reshape(mv,[prod(movieSize(1:3)),movieSize(5)]);
    maxList = max(reshMovie,[],1);
    meanMaxList = mean(maxList);
    stdMaxList = std(maxList);
    %exclude everything that is 4 std above mean
    tooMuchIdx = find(maxList>meanMaxList+4*stdMaxList);
    maxList(tooMuchIdx)=[];
    maxMovie = max(maxList);
    mv = mv/maxMovie;
    
%     %check if already an open panel exists
%     imgFigH = findall(0,'Type','figure','Tag','LabelPanel');
%     if ~isempty(imgFigH)
%         delete(imgFigH);
%     end;
    
    %make max-Projection
    mIntp = max(mv(:,:,:,1,1),[],3);
    
    %open uiViewPanel & show img
    imgFigH = uiViewPanelShowImg(mIntp,1);
    
    %turn off menu 'Panels'
    uiCh = get(imgFigH,'Children');
    set(uiCh(4),'Visible','off');
    
    %set ID tag and name
    set(imgFigH,'Tag','LabelPanel','NumberTitle','off');
    
    %set position
    if ~isempty(labelPanelPos)
        set(imgFigH,'Position',labelPanelPos);
    else
        labelPanelPos = get(imgFigH,'Position');
    end
    
    imgH = findall(imgFigH,'Type','image');
    set(imgH,'CData',mIntp);
    hold on;
    set(imgFigH,'Name',fname);
    %set closereqfcn
    set(imgFigH,'CloseRequestFcn','labelPanelCloseReq','HandleVisibility','Callback');
    
    
    
    %activate controls &set props
    SetUserData(imgFigH,mv,0);
    stackslideH = findall(contFigH,'Tag','slider4');
    stacktextH = findall(contFigH,'Tag','text2');
    timeslideH = findall(contFigH,'Tag','slider3');
    timetextH = findall(contFigH,'Tag','text1');
    %set stackslider and text values
    set(stackslideH,'Value',1);
    set(stackslideH,'Min',0.999999999);
    set(stackslideH,'Max',size(mv,3));
    set(stackslideH,'SliderStep',[1/size(mv,3) 10/size(mv,3)]);
    set(stacktextH,'String','1');
    
    %set timeslider and text values
    set(timeslideH,'Value',1);
    set(timeslideH,'Min',0.9999999);
    set(timeslideH,'Max',size(mv,5));
    set(timeslideH,'SliderStep',[1/size(mv,5) 10/size(mv,5)]);
    set(timetextH,'String','1');
    
    %save directory from which movie was loaded
    SetUserData(imgFigH,fpath,0,'moviePath');
    
else %stop evaluation
    return
end;



% make old additional windows invisible if there are any open
% and remember positions

%3d-viewer
view3DH = findall(0,'Tag','view3Dgui');
if ishandle(view3DH)
    %remember position
    positions.view3DPos = get(view3DH(1),'Position');
    set(view3DH,'Visible','off');
end
view3DfigH = findall(0,'Tag','view3Dfig');
if ~isempty(view3DfigH)
    %remember position
    positions.view3DfigPos = get(view3DfigH(1),'Position');
    set(view3DfigH,'Visible','off');
end

%xz/yz-viewer
zFigH = findall(0,'Tag','XZYZFigure');
if ishandle(zFigH)
    %remember position
    positions.zFigPos = get(zFigH(1),'Position');
    set(zFigH,'Visible','off');
end

%movieDataFigure
dataFigH = findall(0,'Name','CurrentMovieData');    
if ishandle(dataFigH)
    %remember position
    positions.dataFigPos = get(dataFigH(1),'Position');
    set(dataFigH,'Visible','off');
end

%intFigure
figH = findall(0,'Tag','intFig');
if ishandle(figH)%remember position
    positions.intFigPos = get(figH(1),'Position');
    set(figH,'Visible','off');
end

%disFigure
figH = findall(0,'Tag','disFigH');
if ishandle(figH)
    %remember position
    positions.disFigPos = get(figH(1),'Position');
    set(figH,'Visible','off');
end

%save labelPanelPosition
positions.labelPanelPos = labelPanelPos;

%save positions in labelgui
SetUserData(contFigH,positions,1,'positions');

%-------write new PD and new labelPanelH--------
%get old PD string
contFigHandles = guidata(contFigH);
pdStringC = get(contFigHandles.label_chooseWindow_PD,'String');

%find other occurences of the filename
sameIdx = strmatch(fname,pdStringC);

if ~isempty(sameIdx)
    %adjust filename of current window
    fname = [fname, ' #',num2str(length(sameIdx)+1)];
    set(imgFigH,'Name',fname);
end

%add filename to pdString. If there is nothing loaded yet, the first entry
%will read 'no movie loaded'
if strcmp(pdStringC{1},'no movie loaded')
    pdStringC{1}=fname;
else
    upperString = pdStringC(1:end-2);
    lowerString = pdStringC(end-1:end);
    pdStringC = [upperString;{fname};lowerString];
end

%make current figure current selection
set(contFigHandles.label_chooseWindow_PD,'String',pdStringC,'Value',length(pdStringC)-2);

%add handle to list of handles
labelPanelHandles = GetUserData(contFigH,'labelPanelHandles');
labelPanelHandles = [labelPanelHandles,imgFigH];
SetUserData(contFigH,labelPanelHandles,1);

%and save the current handle as current handle to make things easier
SetUserData(contFigH,imgFigH,1,'currentWindow');