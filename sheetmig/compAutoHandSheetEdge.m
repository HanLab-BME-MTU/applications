function compAutoHandSheetEdge(imageFileListPhase,path_t0_auto,path_t0_hand,toDoListInput,movieFormat)

doSave=0;
padZeros=3;
fps=3;
targetDir=[pwd,filesep,'SheetEdge',filesep];
movieFileName=[pwd,filesep,'mov_compAutoHandSheetEdgeMovie'];
if doSave==1
    if ~isdir(targetDir)
        mkdir(targetDir)
    end
end


if nargin<1 || isempty(imageFileListPhase)
    imageFileListPhase=[pwd,filesep,'Phase'];
    try
        imageFileListPhase=getFileListFromFolder(imageFileListPhase);
    catch exception2
        [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
            'Select First phase Image');
        
        if ~ischar(filename) || ~ischar(pathname)
            return;
        end
        
        imageFileListPhase = getFileStackNames([pathname filesep filename]);
    end
elseif isdir(imageFileListPhase)
    imageFileListPhase=getFileListFromFolder(imageFileListPhase);
else
    isValid = 1;
    for frame = 1:numel(imageFileListPhase)
        isValid = isValid && exist(imageFileListPhase{frame}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end

if nargin<2 || isempty(path_t0_auto)
    load('xDetectEdge.mat')
    sheetBndAuto  = sheetBnD;
    sheetEdgeAuto = sheetEdge;
    sheetMaskAuto = sheetMask;
    toDoListAuto  = toDoList;
    clear sheetMask sheetEdge sheetBnD toDoList; 
end

if nargin<3 || isempty(path_t0_hand)
    load('xDetectEdgeHand.mat')
    sheetBndHand  = sheetBnD;
    sheetEdgeHand = sheetEdge;
    sheetMaskHand = sheetMask;
    toDoListHand  = toDoList;
    clear sheetMask sheetEdge sheetBnD toDoList;
end

maxFrame=min(toDoListHand(end),toDoListAuto(end));
listDifFrames=[];
for frame=1:maxFrame
    checkMat=xor(sheetMaskHand(frame).mat,sheetMaskAuto(frame).mat);
    if sum(checkMat(:))>0
        listDifFrames=[listDifFrames frame];
    end
end

if nargin>3 && ~isempty(toDoListInput)
    listDifFrames=toDoListInput;
end


% Read in the phase contrast images (this might take some time):
k=1;
kmax=numel(listDifFrames);
for frame=listDifFrames
    statusText=['Reading in ',kmax,' images'];
    progressText(k/kmax,statusText)
    I{frame}=double(imread(imageFileListPhase{frame}));
end


% Set the position of the current figure to screen size:
scrsz = get(0,'ScreenSize');
h     = figure();
set(h,'Position',scrsz);
pos   = get(h,'Position');

if nargin < 5 || isempty(movieFormat) || strcmp(movieFormat,'mov') == 1 || strcmp(movieFormat,'MOV')
    movieFileName = [movieFileName,'.mov'];
    MakeQTMovie('start',movieFileName);
    MakeQTMovie('framerate',fps);
    movieFormat ='mov';
elseif isempty(movieFormat) || strcmp(movieFormat,'avi') == 1 || strcmp(movieFormat,'AVI')
    movieFileName = [movieFileName,'.avi'];
    movieFormat ='avi';
end

reshow=1;
firstRun=true;
while reshow==1
    figure(h)
    k=1;
    for frame=listDifFrames
        axes('Position',[0 0 1 1]);
        imagesc(I{frame}), title(['Sheet edge of frame: ',num2str(frame)]);
        hold on
        plot(sheetBndAuto(frame).pos(:,2) ,sheetBndAuto(frame).pos(:,1) ,'r-');
        plot(sheetEdgeAuto(frame).pos(:,2),sheetEdgeAuto(frame).pos(:,1),'r.');
        plot(sheetBndHand(frame).pos(:,2) ,sheetBndHand(frame).pos(:,1) ,'g-');
        plot(sheetEdgeHand(frame).pos(:,2),sheetEdgeHand(frame).pos(:,1),'g.');    
        colormap gray;        
        %xlim(imageRange(2,:));
        %ylim(imageRange(1,:));
        %hold on;
        [rowsI,colsI]=size(I{frame});
        text(rowsI*0.05,colsI*0.05,num2str(frame),'Color','white','FontSize',18);
        % textDeltaCoord = min(diff(imageRange,[],2))/20;
        % text(imageRange(2,1)+textDeltaCoord,imageRange(1,1) + textDeltaCoord,num2str(frame),'Color','white','FontSize',18);
        set(h,'Position',pos);
        hold off
        
        if firstRun
            if strcmp(movieFormat,'mov') == 1 || strcmp(movieFormat,'MOV')
                MakeQTMovie('addfigure');
            else
                M(k) = getframe;
            end
        end
        
        if ~firstRun
            pause(0.1);
        end
        
        if doSave==1
            if ~isdir(targetDir)
                mkdir(targetDir)
            end
            saveas(gcf,[targetDir,'sheetEdge',num2str(frame,['%0.',int2str(padZeros),'d']),'.tiff'],'tiffn');
        end
        k=k+1;
    end
    
    if firstRun
        if strcmp(movieFormat,'mov') == 1 || strcmp(movieFormat,'MOV')
            MakeQTMovie('finish');
        elseif strcmp(movieFormat,'avi') == 1 || strcmp(movieFormat,'AVI')
            movie2avi(M,movieFileName,'fps',fps);
        end
    end
    
    reshow=input('Do you want to reshow the movie 1=Yes, 0=No? [1]: ');
    if isempty(reshow)
        reshow=1;
        firstRun=false;
    end     
end
