function [scores,allScores]=fsmWaveAnalysis(analysis,roi,gridSize,time,freq,minValue,tSampling,moveROI)
% fsmWaveAnalysis performs spectral analysis on the results of FSM package
%
% The scores calculated by the kinetic analysis of FSM are assigned to regions of interest
% (a set of subwindows, a polygon, or a rectangle) and the scores over time are analyzed
% for characteristic frequencies.
%
% SYNOPSIS      [scores,allScores]=fsmWaveAnalysis(roi,gridSize,time,minValue,tSampling,moveROI)
%
% REMARK        Run fsmWaveAnalysis through the graphical user interface fsmCenter 
%               ------------------------------------------------------------------
%
% INPUT         analysis    : [ 0|1 0|1 ], turns on and off the first and second analyses,
%                             respectively
%               roi         : one of 'g', 'r', 'p' - define the region of interest for the
%                                                    spectral analysis
%                                    'g' : grid - the whole image is divided into windows 
%                                          of a specified size
%                                    'p' : polygon - the user must draw a polygon
%                                    'r' : rectangle - the user must draw a rectangle
%               gridSize    : [y x] specifies the window dimensions in pixels in the case
%                             that roi is = 'g'
%               time        : number of frames for low-pass filtering in the time domain
%               freq        : number of frames for low-pass filtering in the frequency domain
%               minValue    : minimal power (normalized, pass a value between 0 and 1) for
%                             a peak to be labeled with its period in the power spectrum
%               tSampling   : sampling time - needed to calculate the period of the signal
%                             knowing its frequency
%               moveROI     : [ 0 | 1 ]; if 1, the roi is moved from frame to frame along
%                             the vector field - if loads the saved vector field from a
%                             previous run of specTackle.
%
% OUTPUT        scores      : an [m x t] matrix containing all scores. scores(i,j) contains 
%                             the score for window i at time point j (for roi='g'). If roi 
%                             is 'p' or 'r' i is always = 1.
%
% DEPENDENCES   fsmWaveAnalysis uses { }
%               fsmMain is used by { fsmCenter }
%
% Aaron Ponti, May 15th, 2003

if nargin~=8
    error('Seven input parameter expected');
end

if analysis==[0 0]
    disp('Nothing to calculate.');
    return
end

% Initialize output parameters
scores=0;
allScores=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SELECT ROIs AND KINSCORE FILES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select image
[fName,dirName] = uigetfile(...
    {'*.tif;*.tiff;*.jpg;*.jpeg','Image Files (*.tif,*.tiff,*.jpg,*.jpeg)';
    '*.tif','TIF files (*.tif)'
    '*.tiff','TIFF files (*.tiff)'
    '*.jpg;','JPG files (*.jpg)'
    '*.jpeg;','JPEG files (*.jpeg)'
    '*.*','All Files (*.*)'},...
    'Select image');
if(isa(fName,'char') & isa(dirName,'char'))
   img=nrm(imread([dirName,fName]),1);
   % Store image size
   imgSize=size(img);
else
   return 
end

% Initialize vectors to store polygon vertex coordinates
px=[];py=[];

% Select region of interest
if roi=='p'
    
    % Display image
    h=figure; imshow(img,[]);
    set(h,'Name','Please draw region of interest (polygon)');
    set(h,'NumberTitle','off');


    % Select polygon
    try
        [bw,px,py]=roipoly;
    catch
        close;
        uiwait(msgbox('No polygon was selected. Quitting.','Error','modal'));
        return
    end
        
    % Reformat px and py 
    px=px';
    py=py';
    
    % Close figure
    close(h);

    % Set number of rois
    windows=1;
    
elseif roi=='g'
    
    % Select windows 
    [Gm,Wm]=selectRegion(img,gridSize,gridSize); %round(1.5*gridSize)
    
    % Set number of rois
    windows=size(Wm,1);
    
    % Rearrange Wm as polygons for inpolygon
    for i=1:windows
        py(i,1:5)=[Wm(i,1) Wm(i,1) Wm(i,3) Wm(i,3) Wm(i,1)];
        px(i,1:5)=[Wm(i,2) Wm(i,4) Wm(i,4) Wm(i,2) Wm(i,2)];
    end 
    
elseif roi=='r'
    
    % Display image
    h=figure; imshow(img,[]);
    set(h,'Name','Please draw region of interest (rectangle)');
    set(h,'NumberTitle','off');

    % Select rectangular roi
    try
        rect=getrect(h);
    catch
        close;
        uiwait(msgbox('No polygon was selected. Quitting.','Error','modal'));
        return
    end
    
    px=[rect(1) rect(1)+rect(3) rect(1)+rect(3) rect(1) rect(1)];
    py=[rect(2) rect(2) rect(2)+rect(4) rect(2)+rect(4) rect(2)];
    
    % Close figure
    close(h);

    % Set number of rois
    windows=1;
else
end

% Check whether the selection was successful
if isempty(py)
    disp('No windows selected for the analysis. Quitting.');
    return
end

% Number of vertices of the polygon(s)
len=size(py,2)-1;

% Select corresponding kinScore###.mat
[fName,dirName] = uigetfile(...
    {'*.mat;','Matlab workspaces (*.mat)';
    '*.*','All Files (*.*)'},...
    'Select first kinScore matrix');
if ~(isa(fName,'char') & isa(dirName,'char'))
   return 
end

% Recover all file names from the stack
outFileList=getFileStackNames([dirName,fName]);

% Number of images
n=length(outFileList);
    
% String format for extension
[path,body,no,ext]=getFilenameBody([dirName,fName]);
s=length(no);
strg=sprintf('%%.%dd',s);
    
% Store numeric indices
firstIndex=str2num(no);
indices=[firstIndex:firstIndex+n-1];

% User input needed?
user=0;

% If needed, prepare the series of vector field names
if moveROI==1

    % Create file name
    fileSepPos=findstr(filesep,dirName);
    if length(fileSepPos)>1 & dirName(end)==filesep
        fileSepPos=fileSepPos(end-1); 
        newDirName=dirName(1:fileSepPos);
        indxStr=sprintf(strg,indices(1));
        firstVectorFieldName=[newDirName,'vectors',filesep,'vectors',indxStr,'.mat'];
        % Set flag
        user=0;
    else
        % Set flag
        user=1;
    end
    
    % Try to load it
    if exist(firstVectorFieldName)==2
        
        % Set flag
        user=0;
        
    else
        
        % Set flag
       user=1;
        
    end
    
    if user==0
        
        % Recover all file names from the stack
        vectorOutFileList=getFileStackNames(firstVectorFieldName);
        
    else
        
        % Loading failed, the user must select corresponding vectors###.mat
        [newFName,newDirName] = uigetfile(...
            {'*.mat;','Matlab workspaces (*.mat)';
            '*.*','All Files (*.*)'},...
            'Select first kinScore matrix');
        if ~(isa(newFName,'char') & isa(newDirName,'char'))
            return 
        end
        % Recover all file names from the stack
        vectorOutFileList=getFileStackNames([newDirName,newFName]);
        
    end
    
    % Calculate the center(s) of mass of the polygon(s)
    yCM=sum(py(:,1:len),2)./len;
    xCM=sum(px(:,1:len),2)./len;
        
    % Store the initial position of the polygon(s)
    py0=py;
    px0=px;
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ORDER SCORES TO WINDOWS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize output vector
scores=zeros(windows,n);

% Initialize waitbar
h=waitbar(0,'Please wait while assigning scores...');

% If a polygon leaves the image, mark it to be removed - Initialize
counter=0;
toBeRem=[];

% Start collecting data
for c1=1:n
    
    % If needed, load current vectors
    if moveROI==1 & c1>1;
        
        % Load vectors from the previous frame
        load(char(vectorOutFileList(c1-1)));
        dV=[vectors(:,3)-vectors(:,1) vectors(:,4)-vectors(:,2)];
        
        % Update the center(s) of mass of the polygon(s)
        yCM=sum(py(:,1:len),2)./len;
        xCM=sum(px(:,1:len),2)./len;
        
        % Find the vector closest to the center of mass of the polygon
        D=createDistanceMatrix([yCM xCM],vectors(:,1:2));
        
        % Move all polygons
        for c2=1:size(D,1)
            vIndx=find(D(c2,:)==min(D(c2,:)));
            if isempty(vIndx)
                error('Could not find a valid vector to move the ROI');
            end
            if length(vIndx)>1
                vIndx=vIndx(1);
            end
            
            % Move the current polygon
            py(c2,:)=py(c2,:)+dV(vIndx,1);
            px(c2,:)=px(c2,:)+dV(vIndx,2);
            
        end
    
    end
    
    % Load current kinScore
    load(char(outFileList(c1)));
    indxStr=sprintf(strg,indices(c1));
    eval(['kinScore=kinScore',indxStr,';clear kinScore',indxStr,';']);
        
    % Extract score coordinates
    indx=find(kinScore(:,2)~=0);
    
    if ~isempty(indx)
        
        % The score coordinates are stored into vectors
        y=kinScore(indx,2);
        x=kinScore(indx,3);
        
        for c3=1:windows % All windows
            
            if moveROI==1
                % If a polygon leaves the image, mark it to be removed
                if (any(py(c3,:)<1 | py(c3,:)>imgSize(1)) | any(px(c3,:)<1 | px(c3,:)>imgSize(2)))    
                    counter=counter+1;
                    toBeRem(counter)=c3;
                end
            end
            
            % Extract those score coordinates that belong to the selected polygon
            scoreIndx=find(inpolygon(y,x,py(c3,:),px(c3,:)));
            extrKinScore=kinScore(scoreIndx,:);
            
            % Calculated score for the region
            totScore=sum(extrKinScore(:,4));
            
            % Store calculated score
            scores(c3,c1)=totScore;
            
        end
        
    end
    
    % Update waitbar
    waitbar(c1/n,h);
    
end

% Close waitbar
close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% REARRANGE SCORES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove the first and last time points from scores if analysis was started from image 1
% (those frames have no score)
toRem=0;
if ~any(scores(:,1)~=0)
    scores(:,1)=[];
    toRem=toRem+1;
end
if ~any(scores(:,end)~=0)
    scores(:,end)=[];
    toRem=toRem+1;
end

% Update number of time points
last=n-toRem;

if moveROI==1
    % If some of the polygons left the image, they have to be removed
    toBeRem=unique(toBeRem);
    scores(toBeRem,:)=[];
    py0(toBeRem,:)=[]; px0(toBeRem,:)=[];
    py(toBeRem,:)=[]; px(toBeRem,:)=[];
    yCM(toBeRem)=[]; xCM(toBeRem)=[];

    % Correct the number of windows
    windows=windows-length(toBeRem);

    if isempty(scores)
        disp('No valid windows.');
        return
    end
end

% Concatenate scores
allScores=zeros(1,prod(size(scores)));
for i=1:windows
    allScores(1,(i-1)*last+1:i*last)=scores(i,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPECTRAL ANALYSIS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if analysis(1)==1
    
    % Calculate and plot power spectrum
    [steps,nPower]=plotPowerSpectrum(allScores,time,freq,0);
    labelPowerSpectrum(steps,nPower,minValue,tSampling);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MEAN SCORE PER WINDOWS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if analysis(2)==1
    
    % Plot score evolution only for roi = rectangle or polygon
    if roi=='r' | roi=='p'
        figure;
        plot(scores,'k:','LineWidth',2);
        hold on
        t=[-1:1];g=gauss1d(t,1);
        filtScores=conv(scores,g);
        filtScores=filtScores(3:end-2);
        plot(filtScores,'k-','LineWidth',3);
        title('Score evolution')
        hold off;
    end
    
    % Plot mean score per window
    mWS=mean(scores,2);
    
    % Span
    upperB=max(abs([max(mWS) min(mWS)]));
    span=upperB/15;
    
    % Divide into classes
    mWS=round(mWS/span);
    
    % Convert to pointers to cmap
    mWS(find(mWS<-15))=-15;
    mWS(find(mWS>15))=15;
    mWS=mWS+16;
    
    % Prepare colormap
    cmap=colormap(cMFsm);
    
    % Display image
    h=figure; imshow(img,[]);
    hold on;
    
    if moveROI==1
        % Calculate initial center(s) of mass
        yCM0=sum(py0(:,1:len),2)./len;
        xCM0=sum(px0(:,1:len),2)./len;
    end
    
    for i=1:windows
        
        if moveROI==1
            
            % Draw displacement from the original polygon position
            plot([xCM0(i) xCM(i)],[yCM0(i) yCM(i)],'.-','MarkerSize',24,'LineWidth',3,'Color',cmap(mWS(i),:));

        end
        
        % Display color-coded windows
        fill(px(i,:),py(i,:),cmap(mWS(i),:),'FaceAlpha',0.5);
        
    end
    
    % Put colorbar in another window
    figure;
    scale=repmat([-upperB:upperB/250:upperB],50,1);
    imshow(scale,[]);colormap(cmap);
    title('Colormap: blue: depoly -> red: poly');
    
    hold on
    
end
    