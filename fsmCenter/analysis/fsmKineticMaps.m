function [polyMap,depolyMap,kinMap2C,outputdir]=fsmKineticMaps(projDir,n,sigma)
% fsmKineticMaps creates maps of polymerization, depolymerization and net assembly rate
%
% SYNOPSIS      [polyMap,depolyMap,kinMap2C]=fsmKineticMaps(projDir,n,sigma)
%
% INPUT         projDir        : string pointing to the project directory (where the fsmParam.mat file is 
%                                located). Pass projDir=[] to manually select fsmParam.mat via a dialog.
%               n              : [k m] = [ number of frames to be considered from the whole stack , number of frames for time integration ]
%                                
%                                k = n(1) is the number of frames (actually their corresponding kinScore###.mat files) from the total stack to 
%                                         be considered. If k is 0, the all the frames in the stack are considered. If k=-1, the user will be 
%                                         asked to pick the last kinScore###.mat file at runtime through a dialog.
%                                m = n(2) is the number of frames for temporal averaging of the kinetic scores. If m is 0, then k frames are
%                                         considered. If k is zero as well, the kinetic scores are averaged over the whole stack.
%               sigma (pixels) : sigma for the low-pass filtering of the maps
%                                (default - sigma=5)
%
% OUTPUT        polyMap        : 2D polymerization map integrated over n frames
%               depolyMap      : 2D depolymerization map integrated over n frames
%               netMap         : dual color 2D kinetic map integrated over n frames
%                                red channel for polymerization; green channel for depolymerization
%                                The netMap maps are always normalized with respect to the highest 
%                                score (positive or negative); polyMap and depolyMap are not stretched
%               outputdir      : If n=k, (k any number), many turnover maps will be created in the user-selected range of "kineticScore###.mat" files. 
%                                Only the last polyMap, depolyMap, and netMap will be returned, all other maps will be saved to disk: 'outputdir' is 
%                                the directory where turnover maps are saved (it will be selected by the user at runtime via a dialog).
%                                In contrast, if n=-1,  polyMap, depolyMap, and netMap will be averaged over the whole range, and returned as outputs of the
%                                function. No maps will be saved, and outputdir will be empty.
%
% DEPENDENCES   fsmKineticMaps uses { }
%               fsmKineticMaps is used by {  }
%
% Aaron Ponti, September 2th, 2003

if nargin<2 | nargin>3
    error('Two or three input parameters expected');
end

if nargin==2
    sigma=5;
end

% Initialize outputs
polyMap=[];depolyMap=[];kinMap2C=[];outputdir=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LOOK FOR AND CHECK fsmParam.mat
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fsmParam,status]=fsmPostGetFsmParam(projDir);
if status==0
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHECK THAT THE KINETIC MODULE IS UP-TO-DATE AND GET THE LIST OF kinScore###.mat FILES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(fsmParam.kin,'uptodate')

    % Check that the tracking module is up-to-date
    if fsmParam.kin.uptodate==0
        
        % The tracking module is not up-to-date. Inform the user and return    
        errordlg('The ''kinScore###.mat'' files are not up-to-date. Please re-run the Kinetic module in SpeckTackle.','Error','modal');
        return

    else
        uptodate=1;    
    end
    
else
    
    % Old version of fsmParam.mat. Inform the user that he/she has to make sure that everything is up-to-date
    uiwait(msgbox('Since ''fsmParam.mat'' has been created by an old version of SpeckTackle, I cannot make sure that the ''kinScore###.mat'' files are up-to-date. Continue at your own risk.','Warning'));
    uptodate=-1;
    
end

% Get list of kinScore files
[kinScoreList,success]=fsmPostGetSubProjFileList([projDir,filesep,'kinScore'],'kinScore');
len=length(kinScoreList);

% String format for extension
strg=fsmParam.specific.formString;

% Image size
imgSize=fsmParam.specific.imgSize;

% User input asked
if n(1)==-1

    % Select frames to be processed
    [uFirst,uLast]=fsmTrackSelectFramesGUI(1,len,0,'Select frames to be processed:');

    if uFirst==-1
        return % The user closed the dialog
    end
    
    % Keep only the file names in the user-selected range
    kinScoreList=kinScoreList(uFirst:uLast);
    len=length(kinScoreList); % Update
    
end

% If the user picked the last frame at runtime, set n(1)=length(outFileList)
if n(1)==-1
    n(1)=len;
else
    kinScoreList=kinScoreList(1:n(1));
    len=length(kinScoreList);
end

% If the number of frames for averaging is not set, set it equal to the total number of frames
if n(2)==0 | n(2)>n(1)
    n(2)=n(1);
end

% Check whether files have to be saved to disk
if n(2)==n(1)
    SAVEFILE=0;
else
    SAVEFILE=1;
end

% Select output dir
if SAVEFILE==1
    outputdir=uigetdir('','Select directory to save turnover maps to.');
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
else
    outputdir=[];
end

% Create vector of indices for file names
[path,body,indxStart,ext]=getFilenameBody(char(kinScoreList(1)));
[path,body,indxEnd,ext]=getFilenameBody(char(kinScoreList(end)));
indices=[str2num(indxStart):str2num(indxEnd)-n(2)+1]+fix(n(2)/2);

% Number of images
nImg=length(kinScoreList)-(n(2)-1);

if nImg>1
    
    % Initialize waitbar
    h=waitbar(0,'Creating kinetic maps...');

end

% Create map
for i=1:nImg
    
    % Initialize emtpy maps
    polyMap=zeros(imgSize);
    depolyMap=polyMap;

    for j=1:n(2)
        
        % Current kinScore index
        currentIndx=i+j-1;
        
        % Load kinScore
        load(char(kinScoreList(currentIndx)));
        
        % Read index
        [path,body,indxStr,ext]=getFilenameBody(char(kinScoreList(currentIndx)));
        
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
    
    % Average over n frames
    polyMap=polyMap./n(2);
    depolyMap=depolyMap./n(2);
    
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
    
    if SAVEFILE==1
        
        % Save image and workspaces
        indxStr=sprintf(strg,indices(i));
        fname=[outputdir,filesep,'tif',filesep,'kinMap2C_',indxStr,'.tif'];
        imwrite(kinMap2C,fname,'tif','Compression','none');
        eval(['save ',outputdir,filesep,'mat',filesep,'kinMap2C_',indxStr,'.mat kinMap2C;']);
        eval(['save ',outputdir,filesep,'mat',filesep,'polyMap_',indxStr,'.mat polyMap;']);
        eval(['save ',outputdir,filesep,'mat',filesep,'depolyMap_',indxStr,'.mat depolyMap;']);

    end
    
    if nImg>1

        % Update waitbar
        waitbar(i/nImg,h);
        
    end
    
end

if nImg>1

    % Close waitbar
    close(h)

end
