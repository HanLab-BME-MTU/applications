function fsmSpeedMaps(gridSize,n,d0_init,loadMPM,loadPolygon,sampling,pixelSize,resizeImg,maxSpeed)
% fsmSpeedMaps creates speed maps from the flow maps returned by the SpeckTackle package
%
% fsmSpeedMaps goes through the whole M (or Md) stack of s matrices (each matrix corresponds to the
%    matches returned by the tracker for frames 2 consecutive frames) and creates t<=s speed maps
%    each of which is the average over n frames [j-n/2:j+n/2] around frame j, j in [fix(n/2)+1:s-fix(n/2)]
%   
% SYNOPSIS      fsmSpeedMaps(gridSize,n,d0_init,loadMPM,loadPolygon,sampling,pixelSize,resizeImg,maxSpeed)
%
% INPUT         gridSize    : [y x], defines the grid on which the velocities are calculated   
%               n           : number of frames for temporal averaging.
%                             If n==0, the whole stack is averaged into one speed map.
%               d0_init     : correlation length (pixels)
%               loadMPM     : [ 0 | 1 ] - if 1 loads raw vectors from mpm.mat; if 0 loads
%                                         already filtered vectors (by a previous run of 
%                                         fsmSpeedMaps) from Md_##-##.mat
%               loadPolygon : loads a polygon with which to cut the velocity map at the cell edge
%               sampling    : movie sampling size (s)
%               pixelSize   : pixel size in the image domain (nm)
%               resizeImg   : [ 0 | 1 ] resize the velocity map to original image dimensions
%               maxSpeed    : [ 0 | n ] maximum expected speed (to set the same color scaling for all frames)
%                             Set it to 0 to turn off rescaling (the function will set this value to 110% 
%                             of the maximum velocity from frame 1) or to any velocity n in nm/min.
%
% OUTPUT        velocityMaps saved to disk as 
%
% DEPENDENCES   fsmSpeedMaps uses { }
%               fsmSpeedMaps is used by { }
%
% Aaron Ponti, January 16th, 2004

% Check input parameters

global uFirst uLast

% Load first image
[fName,dirName] = uigetfile(...
    {'*.tif;*.tiff;*.jpg;*.jpeg','Image Files (*.tif,*.tiff,*.jpg,*.jpeg)';
    '*.tif','TIF files (*.tif)'
    '*.tiff','TIFF files (*.tiff)'
    '*.jpg;','JPG files (*.jpg)'
    '*.jpeg;','JPEG files (*.jpeg)'
    '*.*','All Files (*.*)'},...
    'Select one image from the stack');
if(isa(fName,'char') & isa(dirName,'char'))
    img=nrm(imread([dirName,fName]),1);
    % Store image size
    imgSize=size(img);
else
    return 
end

% String format for extension
[path,outputFileName,no,ext]=getFilenameBody([dirName,fName]);
s=length(no);
strg=sprintf('%%.%dd',s);

if loadMPM==1
    
    % Load MPM.mat
    [fName,dirName] = uigetfile(...
        {'MPM.mat;','MPM.mat';
        '*.*','All Files (*.*)'},...
        'Select MPM.mat');
    if ~(isa(fName,'char') & isa(dirName,'char'))
        return 
    end
    load([dirName,fName]);
    clear MPM;
    
    % Toggle interpolation
    interp=1;
    
else
    
    % Load Md.mat
    [fName,dirName] = uigetfile(...
        {'Md.mat;','Md.mat';
        '*.*','All Files (*.*)'},...
        'Select Md.mat');
    if ~(isa(fName,'char') & isa(dirName,'char'))
        return 
    end
    load([dirName,fName]);
    
    % No interpolation
    interp=0;

end

if loadPolygon==1
    % Load corresponding polygon
    [fName,dirName] = uigetfile(...
        {'*.mat;','Matlab workspaces (*.mat)';
        '*.*','All Files (*.*)'},...
        'Select polygon.mat - Click cancel to skip');
    if ~(isa(fName,'char') & isa(dirName,'char'))
        polygon=[]; 
    else
        load([dirName,fName]);
    end
else
    polygon=[];
end

% Select first and last frame
frames=size(M,3);
% The minimum number of frames to be selected depends on the value assigned to n by the user
if n==0
    minN=1; % This causes the all stack to be averaged into one speed map
else
    minN=n-1;
end

guiH=fsmTrackSelectFramesGUI; ch=get(guiH,'Children');
set(findobj('Tag','pushOkay'),'UserData',minN); % At least n-1 frames must be considered
title='Select frame pairs to be processed:';
set(findobj('Tag','editFirstFrame'),'String',num2str(1));
set(findobj('Tag','editLastFrame'),'String',num2str(frames));
set(findobj('Tag','SelectFramesGUI'),'Name',title);
sSteps=[1/(frames-1) 1/(frames-1)];
set(findobj('Tag','sliderFirstFrame'),'SliderStep',sSteps,'Max',frames,'Min',1,'Value',1);
set(findobj('Tag','sliderLastFrame'),'SliderStep',sSteps,'Max',frames,'Min',1,'Value',frames);
waitfor(guiH); % The function waits for the dialog to close (and to return values for uFirst and uLast)

if uFirst==-1
    return % The user closed the dialog
end

% Update n depending on the selection
if n==0
    n=uLast;
end

% Interpolation grid
G=framework(imgSize,gridSize);

% Select output dir
outputdir=uigetdir('','Select directory for output');

% Start
if interp==1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % PHASE 1: INTERPOLATION
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    h=waitbar(0,'Interpolating...');
    
    % Interpolate vector fields through all time points
    current=0; 
    steps=1+uLast-uFirst;
    for c1=uFirst:uLast
        
        % Update frame counter
        current=current+1;
        
        % Extract vectors
        Mp=M(:,:,c1);
        Mv=Mp(find(Mp(:,1)~=0 & Mp(:,3)~=0),:);
        
        % Interpolate onto a grid
        [divM,d0]=vectorFieldDiv(Mv,G,d0_init,polygon);
        d0=updateD0FromDiv(divM,d0_init,1,size(Mv,1),size(G,1));
        Md(:,:,current)=vectorFieldInterp(Mv,G,d0,polygon);
        
        % Updating waitbar
        waitbar(current/steps,h);
        
    end
    
    % Save interpolated vector fields (Md) to disk
    eval(['save ',outputdir,filesep,'Md_',sprintf(strg,uFirst),'-',sprintf(strg,uLast),'_d0=',num2str(d0_init),'.mat Md;']);
    
    % Close waitbar
    close(h);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PHASE 2: TIME AVERAGING AND FIGURE CREATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create vector of indices for file names
indices=[uFirst:uLast-n+1]+fix(n/2);

h=waitbar(0,'Creating speed maps...');

% Calculate average vector fields
steps=size(Md,3)-n+1;
for c2=1:steps
    
    % Extract mean velocities over n frames
    vectors=Md(:,3:4,c2:c2+n-1)-Md(:,1:2,c2:c2+n-1);
    meanVectors=mean(vectors,3);
    
    % Construct average vector matrix
    Mav=zeros(size(Md,1),size(Md,2));
    Mav(:,1:2)=Md(:,1:2,1);
    Mav(:,3:4)=Mav(:,1:2)+meanVectors;
    
    % Calculate velocities
    velocities=zeros(size(Mav,1),3);
    velocities(:,1:2)=Mav(:,1:2);
    velocities(:,3)=sqrt(meanVectors(:,1).^2+meanVectors(:,2).^2);
    
    % Reshape
    velocityMap=reshape(velocities(:,3),length(unique(G(:,2))),length(unique(G(:,1))))';
    
    % Transform to nm/min
    velocityMap=velocityMap*(60/sampling)*pixelSize;
    
    % If needed, resize
    if resizeImg==1
        velocityMap=imresize(velocityMap,imgSize,'bicubic');
    end
    
    % Display
    figure
    imshow(velocityMap,[]);
    colormap('jet')
    colorbar;
    
    % Set color range between 0 and MAXSPEED and update colorbar
    if c2==1 & maxSpeed==0
        maxSpeed=max(velocityMap(:))*1.1; % Give a 10% space
    end
    set(gca,'CLim',[0 maxSpeed]);
    colorbar;
    
    % Save image
    indxStr=sprintf(strg,indices(c2));
    fname=[outputdir,filesep,outputFileName,indxStr,'.tif'];
    saveas(gcf,fname,'tif');
    
    % Save velocityMap to disk as well
    eval(['save ',outputdir,filesep,outputFileName,'velocityMap_d0=',num2str(d0_init),'_',indxStr,'.mat velocityMap']);
    
    close(gcf);
    
    % Updating waitbar
    waitbar(c2/steps,h);
    
end

close(h);
