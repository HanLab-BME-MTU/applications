function outputdir=fsmSpeedMaps(gridSize,n,d0_init,loadMPM,sampling,pixelSize,overlayVect,userROIbw,maxSpeed,segment,bitDepth)
% fsmSpeedMaps creates speed maps from the flow maps returned by the SpeckTackle package
%
% fsmSpeedMaps goes through the whole M (or Md) stack of s matrices (each matrix corresponds to the
%    matches returned by the tracker for frames 2 consecutive frames) and creates t<=s speed maps
%    each of which is the average over n frames [j-n/2:j+n/2] around frame j, j in [fix(n/2)+1:s-fix(n/2)]
%   
% SYNOPSIS      outputdir=fsmSpeedMaps(gridSize,n,d0_init,loadMPM,sampling,pixelSize,overlayVect,userROIbw,maxSpeed)
%
% INPUT         gridSize    : [y x], defines the grid on which the velocities are calculated   
%               n           : number of frames for temporal averaging.
%                             If n==0, the whole stack is averaged into one speed map.
%               d0_init     : correlation length (pixels)
%               loadMPM     : [ 0 | 1 ] - if 1 loads raw vectors from mpm.mat; if 0 loads
%                                         already filtered vectors (by a previous run of 
%                                         fsmSpeedMaps) from Md_##-##.mat
%               sampling    : movie sampling size (s)
%               pixelSize   : pixel size in the image domain (nm)
%               overlayVect : [ 0 | 1 ] overlays vector field to speed map
%               userROIbw   : black and white mask (0,1)-matrix used to mask parts of the image
%               maxSpeed    : [ 0 | n ] maximum expected speed (to set the same color scaling for all frames)
%                             Set it to 0 to turn off rescaling (the function will set this value to 110% 
%                             of the maximum velocity from frame 1) or to any velocity n in nm/min.
%               segment     : [ 0 | 1 ] turns on and off automatic segmentation
%               bitDpeth    : image bit depth, used for automatic segmentation
%
% OUTPUT        outputdir   : directory where the speed maps are saved to.
%               speedMaps saved to disk as .tif, .eps, .mat 
%
%
% DEPENDENCES   fsmSpeedMaps uses { }
%               fsmSpeedMaps is used by { }
%
% Aaron Ponti, January 16th, 2004

global uFirst uLast

% Check input parameters
if nargin~=11
    error('11 input parameters expected');
end

% Initialize output
outputdir=[];

% Load first image
[fName,dirName] = uigetfile(...
    {'*.tif;*.tiff;*.jpg;*.jpeg','Image Files (*.tif,*.tiff,*.jpg,*.jpeg)';
    '*.tif','TIF files (*.tif)'
    '*.tiff','TIFF files (*.tiff)'
    '*.jpg;','JPG files (*.jpg)'
    '*.jpeg;','JPEG files (*.jpeg)'
    '*.*','All Files (*.*)'},...
    'Select first image');
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

% Read image file list
outFileList=getFileStackNames([dirName,fName]);

if loadMPM==1
    
    % Load mpm.mat
    [fName,dirName] = uigetfile(...
        {'mpm.mat;','mpm.mat';
        '*.*','All Files (*.*)'},...
        'Select mpm.mat');
    if ~(isa(fName,'char') & isa(dirName,'char'))
        return 
    end
    load([dirName,fName]);
    if ~exist('M')
        % The user didn't pick a valid mpm.mat file
        uiwait(msgbox('Invalid mpm.mat file chosen.','error','modal'));
        return
    end
    clear MPM;
    
    % Toggle interpolation
    interp=1;
    
    % Number of frames
    frames=size(M,3);

else
    
    % Load Md.mat
    [fName,dirName] = uigetfile(...
        {'*.mat;','*.mat';
        '*.*','All Files (*.*)'},...
        'Select Md###-###.mat');
    if ~(isa(fName,'char') & isa(dirName,'char'))
        return 
    end
    load([dirName,fName]);
    if ~exist('Md')
        % The user didn't pick a valid Md.mat file
        uiwait(msgbox('Invalid Md###-###.mat file chosen.','error','modal'));
        return
    end

    % No interpolation
    interp=0;

    % Number of frames
    frames=size(Md,3);

end

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

if n>uLast
    n=uLast;
end

if loadMPM==1
    % Interpolation grid
    G=framework(imgSize,gridSize);
end

% Select output dir
outputdir=uigetdir('','Select directory to save speed maps to.');
if outputdir==0 % The user clicked on cancel
    disp('Aborted by the user.');
    outputdir=[];
    return
end

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
        [divM,d0]=vectorFieldDiv(Mv,G,d0_init,[]);
        d0=updateD0FromDiv(divM,d0_init,1,size(Mv,1),size(G,1));
        Md(:,:,current)=vectorFieldInterp(Mv,G,d0,[]);
        
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

% Create subdirectories if needed
if exist([outputdir,filesep,'tif'])~=7
    % Create directory
    success=mkdir(outputdir,'tif');
    if success==0
        error('Could not create subfolder in specified directory');
    end
end
if exist([outputdir,filesep,'mat'])~=7
    % Create directory
    success=mkdir(outputdir,'mat');
    if success==0
        error('Could not create subfolder in specified directory');
    end

end
if exist([outputdir,filesep,'eps'])~=7
    % Create directory
    success=mkdir(outputdir,'eps');
    if success==0
        error('Could not create subfolder in specified directory');
    end
end

% Create vector of indices for file names
[path,body,indxStart,ext]=getFilenameBody(char(outFileList(1)));
[path,body,indxEnd,ext]=getFilenameBody(char(outFileList(end)));
indices=[str2num(indxStart):str2num(indxEnd)-n+1]+fix(n/2);

% Initializing waitbar
h=waitbar(0,'Creating speed maps...');

% Create a full bwMask (in case of failure of edge detection)
bwMask=ones(size(img));

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
    speedMap=reshape(velocities(:,3),length(unique(Md(:,2,1))),length(unique(Md(:,1,1))))';
    
    % Transform to nm/min
    speedMap=speedMap*(60/sampling)*pixelSize;
    
    % If needed, resize
    speedMap=imresize(speedMap,imgSize,'bilinear');
    
    if segment==1
        % Find cell boundaries
        xmax=2^bitDepth-1;
        img=imreadnd2(char(outFileList(c2)),0,xmax);
        try
            [ans,img_edge,bwMask]=imFindCellEdge(img,'',0,'filter_image',1,'img_sigma',1,'bit_depth',xmax);
        catch
            % Uses last one
        end       
    
        % Crop velocity map
        speedMap=speedMap.*bwMask;
        
    end
    
    % If needed apply the userROI as well
    if ~isempty(userROIbw)
        % Check size
        if size(speedMap)==size(userROIbw)
            speedMap=speedMap.*userROIbw;
        end
    end
    
    % Display
    mH=figure;
    imshow(speedMap,[]);
    colormap('jet')

    % Plot vectors if requested by the user
    if overlayVect==1
        
        % Remove vectors outside the mask
        toBeDispl=[];
        counter=0;
        for i=1:size(Mav,1)
            if bwMask(Mav(i,1),Mav(i,2))==1
                counter=counter+1;
                toBeDispl(counter)=i;
            end
        end
        MavDisp=Mav(toBeDispl,:);
        
        % Scale factor
        if c2==1
            scaleFactor=fix(min(gridSize)/mean(sqrt((Mav(:,3)-Mav(:,1)).^2+(Mav(:,4)-Mav(:,2)).^2)));
        end
        
        % Plot vectors on top of speed map
        MavDisp(:,3:4)=[MavDisp(:,1:2)+scaleFactor*(MavDisp(:,3:4)-MavDisp(:,1:2))];
        hold on;
        % Check for MATLAB version - quiver 7.0 is no longer compatible
        v=ver('MATLAB');
        pointPos=findstr(v.Version,'.');
        if ~isempty(pointPos)
            v.Version=v.Version(1:pointPos(1)+1);
        end
        if str2num(v.Version)<7
            qH=quiver(MavDisp(:,2),MavDisp(:,1),MavDisp(:,4)-MavDisp(:,2),MavDisp(:,3)-MavDisp(:,1),0);
        else
            qH=quiver('v6',MavDisp(:,2),MavDisp(:,1),MavDisp(:,4)-MavDisp(:,2),MavDisp(:,3)-MavDisp(:,1),0);
        end
        set(qH(1),'Color',[0 0 0]); % Base
        set(qH(2),'Color',[0 0 0]); % Tip
    else
        scaleFactor=0;
    end

    % Set color range between 0 and MAXSPEED and update colorbar
    if c2==1 & maxSpeed==0
        maxSpeed=max(speedMap(:))*1.1; % Give a 10% space
    end
    set(gca,'CLim',[0 maxSpeed]);
    colorbar;

    % Save image
    indxStr=sprintf(strg,indices(c2));
    if scaleFactor~=0
        fname=[outputdir,filesep,'tif',filesep,'speedMap_d0=',num2str(d0_init),'_scale=',num2str(scaleFactor),'x_frames=',num2str(n),'_',indxStr,'.tif'];
        print(gcf,'-dtiffnocompression',fname);
        fname=[outputdir,filesep,'eps',filesep,'speedMap_d0=',num2str(d0_init),'_scale=',num2str(scaleFactor),'x_frames=',num2str(n),'_',indxStr,'.eps'];
        print(gcf,'-dpsc2',fname);
        % Save speedMap to disk as well
        eval(['save ',outputdir,filesep,'mat',filesep,'speedMap_d0=',num2str(d0_init),'_scale=',num2str(scaleFactor),'x_frames=',num2str(n),'_',indxStr,'.mat speedMap;']);
    else
        fname=[outputdir,filesep,'tif',filesep,'speedMap_d0=',num2str(d0_init),'_frames=',num2str(n),'_',indxStr,'.tif'];
        print(gcf,'-dtiffnocompression',fname);
        fname=[outputdir,filesep,'eps',filesep,'speedMap_d0=',num2str(d0_init),'_frames=',num2str(n),'_',indxStr,'.eps'];
        print(gcf,'-dpsc2',fname);
        % Save speedMap to disk as well
        eval(['save ',outputdir,filesep,'mat',filesep,'speedMap_d0=',num2str(d0_init),'_frames=',num2str(n),'_',indxStr,'.mat speedMap;']);
    end
    


    close(gcf);
    
    % Updating waitbar
    waitbar(c2/steps,h);
    
end

close(h);
