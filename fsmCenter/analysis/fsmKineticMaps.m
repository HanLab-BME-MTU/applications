function [polyMap,depolyMap,kinMap2C]=fsmKineticMaps(firstKinScore,imgSize,n,sigma)
% fsmKineticMaps creates maps of polymerization, depolymerization and net assembly rate
%
% SYNOPSIS      [polyMap,depolyMap,kinMap2C]=fsmKineticMaps(firstKinScore,imgSize,n,sigma)
%
% INPUT         firstKinScore  : string containing the name (with complete path) of the 
%                                first kinScore###.mat file
%                                set firstKinScore=[] to have the function open a dialog
%               imgSize        : size of the analyzed images
%               n              : number of frames for time integration (0 for the entire stack)
%               sigma (pixels) : sigma for the low-pass filtering of the maps
%                                (default - sigma=5)
%
% OUTPUT        polyMap        : 2D polymerization map integrated over n frames
%               depolyMap      : 2D depolymerization map integrated over n frames
%               netMap         : dual color 2D kinetic map integrated over n frames
%                                red channel for polymerization; green channel for depolymerization
%                                The netMap maps are always normalized with respect to the highest 
%                                score (positive or negative); polyMap and depolyMap are not stretched
%                                All maps are written to disk.
%
% DEPENDENCES   fsmKineticMaps uses { }
%               fsmKineticMaps is used by {  }
%
% Aaron Ponti, September 2th, 2003

global uFirst uLast

if nargin<3 | nargin>4
    error('Three or four input parameters expected');
end

if nargin==3
    sigma=5;
end

if isempty(firstKinScore) | exist(firstKinScore)~=2 % Not a file
    
    % Select kinScore###.mat
    [fName,dirName] = uigetfile(...
        {'*.mat;','Matlab workspaces (*.mat)';
        '*.*','All Files (*.*)'},...
        'Select first kinScore matrix');
    if ~(isa(fName,'char') & isa(dirName,'char'))
        polyMap=[];depolyMap=[];img2C=[];
        return 
    end
    
    firstKinScore=[dirName,fName];
end

% String format for extension
[path,outputFileName,no,ext]=getFilenameBody([dirName,fName]);
s=length(no);
strg=sprintf('%%.%dd',s);

% Recover all file names from the stack
outFileList=getFileStackNames(firstKinScore);
len=length(outFileList);

% Select range of frames for which to create maps
guiH=fsmTrackSelectFramesGUI; ch=get(guiH,'Children');
set(findobj('Tag','pushOkay'),'UserData',0); % Minimum range 
title='Select images to be processed:';
set(findobj('Tag','editFirstFrame'),'String',num2str(1));
set(findobj('Tag','editLastFrame'),'String',num2str(len));
set(findobj('Tag','SelectFramesGUI'),'Name',title);
sSteps=[1/(len-1) 1/(len-1)];
set(findobj('Tag','sliderFirstFrame'),'SliderStep',sSteps,'Max',len,'Min',1,'Value',1);
set(findobj('Tag','sliderLastFrame'),'SliderStep',sSteps,'Max',len,'Min',1,'Value',len);
waitfor(guiH); % The function waits for the dialog to close (and to return values for uFirst and uLast)

if uFirst==-1
    return % The user closed the dialog
end

% Keep only the file names in the user-selected range
outFileList=outFileList(uFirst:uLast);
len=length(outFileList); % Update

% Only consider first n frames if n~=0
if n==0
    n=len;
end

% Create vector of indices for file names
[path,body,indxStart,ext]=getFilenameBody(char(outFileList(1)));
[path,body,indxEnd,ext]=getFilenameBody(char(outFileList(end)));
indices=[str2num(indxStart):str2num(indxEnd)-n+1]+fix(n/2);

% Select output dir
outputdir=uigetdir('','Select directory for output');
if outputdir==0 % The user clicked on cancel
    disp('Aborted by the user.');
    return
end

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

% Initialize waitbar
h=waitbar(0,'Creating kinetic maps...');

% Number of images
nImg=length(outFileList)-(n-1);

% Create map
for i=1:nImg
    
    % Initialize emtpy maps
    polyMap=zeros(imgSize);
    depolyMap=polyMap;

    for j=1:n
        
        % Current kinScore index
        currentIndx=i+j-1;
        
        % Load kinScore
        load(char(outFileList(currentIndx)));
        
        % Read index
        [path,body,indxStr,ext]=getFilenameBody(char(outFileList(currentIndx)));
        
        % Copy matrix
        eval(['kinScore=kinScore',indxStr,'; clear kinScore',indxStr,';']);
        
        % Find scores
        indx=find(kinScore(:,2)~=0);
        
        % Copy scores into map
        if ~isempty(indx)
            
            for k=1:length(indx)
                
                % Read score
                score=kinScore(k,4);
                
                switch sign(score)
                    case 1, polyMap(kinScore(k,2),kinScore(k,3))=polyMap(kinScore(k,2),kinScore(k,3))+score;
                    case -1, depolyMap(kinScore(k,2),kinScore(k,3))=depolyMap(kinScore(k,2),kinScore(k,3))+score;
                    otherwise
                        error('A zero score.');
                end
                
            end
            
        end
        
    end
    
    % Low-pass filter
    if sigma~=0
        polyMap=Gauss2D(polyMap,sigma);
        depolyMap=Gauss2D(depolyMap,sigma);
    end
    
    % Create dual-channel image
    mx=max([polyMap(:);abs(depolyMap(:))]);
    
    kinMap2C=zeros([imgSize 3]);
    kinMap2C(:,:,1)=polyMap/mx;
    kinMap2C(:,:,2)=abs(depolyMap)/mx;
    
    % Save image and workspaces
    indxStr=sprintf(strg,indices(i));
    fname=[outputdir,filesep,'tif',filesep,'kinMap2C_',indxStr,'.tif'];
    imwrite(kinMap2C,fname,'tif','Compression','none');
    eval(['save ',outputdir,filesep,'mat',filesep,'kinMap2C_',indxStr,'.mat kinMap2C;']);
    eval(['save ',outputdir,filesep,'mat',filesep,'polyMap_',indxStr,'.mat polyMap;']);
    eval(['save ',outputdir,filesep,'mat',filesep,'depolyMap_',indxStr,'.mat depolyMap;']);

    % Update waitbar
    waitbar(i/nImg,h);

end

% Close waitbar
close(h)

