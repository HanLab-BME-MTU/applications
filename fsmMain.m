function fsmParam=fsmMain(fsmParam)
% fsmMain Main function of the fsm project
%
% SYNOPSIS      fsmMain
%
% REMARK        Run fsmMain through the graphical user interface fsmGuiMain
%               -----------------------------------------------------------
%               
% INPUT         fsmParam    : parameter structure
%                             (type 'help fsmGetparamsDeflt' for more detail)
%
% OUTPUT        fsmParam    : parameter structure completed with experiment specific data 
%
% DEPENDENCES   fsmMain uses { fsmToolsPrepareWorkDir
%                              fsmGetParamDflts
%                              fsmPrepMain
%                              fsmTrackMain
%                              fsmBuildMain
%                              fsmKinMain
%                              fsmDispMain
%                              fsmInfoSpeckleArrayStatistics
%                              }
%               fsmMain is used by { fsmMainRun }
%
% Aaron Ponti, October 2nd, 2002

% TEST

if nargin~=1
    error('One parameter (fsmParam) expected');
end

% Store current directory
olddir=cd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   CREATE WORK DIRECTORY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

status=fsmToolsPrepareWorkDir(fsmParam.main.path);
if status==0
    return;
end

% Make sure to be in the work path
cd(fsmParam.main.path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   START TIMER
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

% Write as well to 'log.txt'
fid=fopen('log.txt','a');

% Could the file be opened successfully?
if fid==-1
	error('Couldn''t open the file.');
end

% Write to file
fprintf(fid,'****************************************************************************************\n');
fprintf(fid,'*                                                                                      *\n');
fprintf(fid,'*     Calculation started: %s                                        *\n',datestr(now,0));
fprintf(fid,'*                                                                                      *\n');
fprintf(fid,'****************************************************************************************\n');

% Close the file
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   PREPROCESSING MODULE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run if enabled
if fsmParam.prep.enable==1
    % Call main function of the module
    [fsmParam,status]=fsmPrepMain(fsmParam);
else
    status=2;
end

% Check the status returned by the module
if status==1
    disp('The PREPROCESSING module completed successfully.');
    
    % Save fsmParam to disk
    eval(['save ',fsmParam.main.path,filesep,'fsmParam.mat fsmParam']);

elseif status==0
    disp('The PREPROCESSING module did not complete successfully. Aborting.');
    return
else
    disp('The PREPROCESSING module was not run.');
end

% If preprocessing is not run, check that the list of images is correct
if status==2
    [fsmParam,status]=fsmCheckFsmParam(fsmParam);
    if status==0
        return;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   TRACKER MODULE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run if enabled
if fsmParam.track.enable==1
    
    if status==2 % The preprocessing module was not run

        % Let the user select the images to be tracked
        [fsmParam,statusSF]=fsmTrackSelectFrames(fsmParam);
        
        if statusSF==0 % Error returned by fsmTrackSelectFrames
            disp('The TRACKING module was interrupted because of an error.');
            return
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   CHECK THAT ALL INPUT PARAMETERS EXIST 
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Check that n (number of images) files cands###.mat exist on disk
    %   It only checks for the existance of the first and the last file
    
    if fsmParam.track.tracker<=2
        if fsmParam.specific.imageNumber<2
            uiwait(msgbox('At least 2 preprocessed images are needed for tracking. Please run the PREPROCESSING module again.','Error','error','modal'));
            return 
        end    
    elseif fsmParam.track.tracker==3
        if fsmParam.specific.imageNumber<3
            uiwait(msgbox('At least 3 preprocessed images are needed for tracking. Please run the PREPROCESSING module again.','Error','error','modal'));
            return
        end
    else
        error('Bad tracker selection in fsmParam.track.tracker');
    end
    
    strg=fsmParam.specific.formString;  % Format string for the correct numbering of the files
    userPath=fsmParam.main.path;
    firstIndex=fsmParam.specific.firstIndex;
    lastIndex=fsmParam.specific.lastIndex;

    % Files to check
    firstFileName=[userPath,filesep,'cands',filesep,'cands',sprintf(strg,firstIndex),'.mat'];
    lastFileName=[userPath,filesep,'cands',filesep,'cands',sprintf(strg,lastIndex),'.mat'];
    if exist(firstFileName)~=2 | exist(lastFileName)~=2 % At least one of the files does not exist
        uiwait(msgbox('Files: ''cands###.mat'' not found. Please run the PREPROCESSING module.','Error','error','modal'));
        return    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Call main function of the module
    [fsmParam,status]=fsmTrackMain(fsmParam);
    
else
    status=2;
end

% Check the status returned by the module
if status==1
    disp('The TRACKING module completed successfully.');

    % Save fsmParam to disk
    eval(['save ',fsmParam.main.path,filesep,'fsmParam.mat fsmParam']);
    
elseif status==0
    disp('The TRACKING module did not complete successfully. Aborting.');
    return
else
    disp('The TRACKING module was not run.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   SPECKLEARRAY BUILDER MODULE
%
%   From the BUILDER module on, at least 4 frames need to have been processed and tracked. 
%   The relaxation of having less than 4 frames processed was introduced to allow the user 
%   to test the outcome of the processing only on 1 frame or to test the tracking only on a 
%   frame pair. For the builder and following modules, though, this is not enough and
%   the software will therefore inform the user and quit.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run if enabled
if fsmParam.build.enable==1

    if size(fsmParam.specific.fileList,1)<4
        uiwait(msgbox('Only the PREPROCESSING and the TRACKER modules are allowed to work with less than 4 frames preprocessed. Please run the PREPROCESSING MODULE again on more frames.','Error','error','modal'));
        disp('Please run the PREPROCESSING module on at least 4 frames. Aborting');
        return 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   CHECK THAT ALL INPUT PARAMETERS EXIST 
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Check that: cands###.mat, gapList###.mat, mpm.mat exist on disk
    
    strg=fsmParam.specific.formString;  % Format string for the correct numbering of the files
    userPath=fsmParam.main.path;
    firstIndex=fsmParam.specific.firstIndex;
    lastIndex=fsmParam.specific.lastIndex;
    
    % Check cands###.mat
    firstFileName=[userPath,filesep,'cands',filesep,'cands',sprintf(strg,firstIndex),'.mat'];
    lastFileName=[userPath,filesep,'cands',filesep,'cands',sprintf(strg,lastIndex),'.mat'];
    if exist(firstFileName)~=2 | exist(lastFileName)~=2 % At least one of the files does not exist
        uiwait(msgbox('Files: ''cands###.mat'' not found. Please run the PREPROCESSING module.','Error','error','modal'));
        return    
    end

    % Check gapList###.mat
    firstFileName=[userPath,filesep,'gapList',filesep,'gapList',sprintf(strg,firstIndex),'.mat'];
    lastFileName=[userPath,filesep,'gapList',filesep,'gapList',sprintf(strg,lastIndex),'.mat'];
    if exist(firstFileName)~=2 | exist(lastFileName)~=2 % At least one of the files does not exist
        uiwait(msgbox('Files: ''gapList###.mat'' not found. Please run the TRACKING module.','Error','error','modal'));
        return    
    end

    % Check mpm.mat
    firstFileName=[userPath,filesep,'mpm.mat'];
    if exist(firstFileName)~=2 % mpm.mat does not exist
        uiwait(msgbox('Files: ''mpm.mat'' not found. Please run the TRACKING module.','Error','error','modal'));
        return    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Call main function of the module
    [fsmParam,status,speckleArray]=fsmBuildMain(fsmParam);
else
    status=2;
end

% Check the status returned by the module
if status==1
    disp('The BUILDER module completed successfully.');
    
    % Save fsmParam to disk
    eval(['save ',fsmParam.main.path,filesep,'fsmParam.mat fsmParam']);

elseif status==0
    disp('The BUILDER module did not complete successfully. Aborting.');
    return
else
    disp('The BUILDER module was not run.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   KINETIC ANALYSIS MODULE
%
%   From the BUILDER module on, at least 4 frames need to have been processed and tracked. 
%   The relaxation of having less than 4 frames processed was introduced to allow the user 
%   to test the outcome of the processing only on 1 frame or to test the tracking only on a 
%   frame pair. For the builder and following modules, though, this is not enough and
%   the software will therefore inform the user and quit.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run if enabled
if fsmParam.kin.enable==1
    
    if size(fsmParam.specific.fileList,1)<4
        uiwait(msgbox('Only the PREPROCESSING and the TRACKER modules are allowed to work with less than 4 frames preprocessed. Please run the PREPROCESSING MODULE again on more frames.','Error','error','modal'));
        disp('Please run the PREPROCESSING module on at least 4 frames. Aborting.');
        return 
    end

    if exist('speckleArray')~=1  
        
        % If the variable speckleArray does not exist, it means that the
        % previous module was not run. So speckleArray.mat has to be loaded 
        % from disk
                               
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %   CHECK THAT ALL INPUT PARAMETERS EXIST 
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Check that: speckleArray.mat exists on disk
        firstFileName=[fsmParam.main.path,filesep,'speckleArray.mat'];
        if exist(firstFileName)~=2  % speckleArray.mat does not exist
            uiwait(msgbox('File ''speckleArray.mat'' not found. Please run the BUILDER module.','Error','error','modal'));
            return    
        else
            load(firstFileName); % Load speckleArray.mat from disk
            
            % Check that the speckleArray is really in the old format
            if length(speckleArray)>1 & length(speckleArray(1).timepoint)==1
                uiwait(msgbox('The loaded speckleArray has the old structure. Please run "Convert specklArray" in fsmCenter and try again.','Error','modal'));
                return
            end

        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Call main function of the module
    [fsmParam,status,speckleArray,SCORE]=fsmKinMain(fsmParam,speckleArray);
    
    % Some statistics
    fsmInfoSpeckleArrayStatistics(speckleArray,1);
    
else
    status=2;
end

% Check the status returned by the module
if status==1
    disp('The KINETIC ANALYSIS module completed successfully.');
    
    % Save fsmParam to disk
    eval(['save ',fsmParam.main.path,filesep,'fsmParam.mat fsmParam']);

elseif status==0
    disp('The KINETIC ANALYSIS module did not complete successfully. Aborting.');
    return
else
    disp('The KINETIC ANALYSIS module was not run.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   RESULT DISPLAY MODULE
%
%   From the BUILDER module on, at least 4 frames need to have been processed and tracked. 
%   The relaxation of having less than 4 frames processed was introduced to allow the user 
%   to test the outcome of the processing only on 1 frame or to test the tracking only on a 
%   frame pair. For the builder and following modules, though, this is not enough and
%   the software will therefore inform the user and quit.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run if enabled
if fsmParam.disp.enable==1
    
    if size(fsmParam.specific.fileList,1)<4
        uiwait(msgbox('Only the PREPROCESSING and the TRACKER modules are allowed to work with less than 4 frames preprocessed. Please run the PREPROCESSING MODULE again on more frames.','Error','error','modal'));
        disp('Please run the PREPROCESSING module on at least 4 frames. Aborting.');
        return 
    end
    
    if exist('SCORE')~=1
        
        % If the variable SCORE does not exist, it means that the
        % previous module was not run. So SCORE.mat has to be loaded 
        % from disk
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %   CHECK THAT ALL INPUT PARAMETERS EXIST 
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Check that: SCORE.mat exists on disk
        firstFileName=[fsmParam.main.path,filesep,'SCORE.mat'];
        if exist(firstFileName)~=2  % speckleArray.mat does not exist
            uiwait(msgbox('File ''SCORE.mat'' not found. Please run the BUILDER module.','Error','error','modal'));
            return    
        else
            load(firstFileName); % Load SCORE.mat from disk
        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Call main function of the module
    [fsmParam,status]=fsmDispMain(fsmParam,SCORE);
else
    status=2;
end

% Check the status returned by the module
if status==1
    disp('The RESULT DISPLAY module completed successfully.');
    
    % Save fsmParam to disk
    eval(['save ',fsmParam.main.path,filesep,'fsmParam.mat fsmParam']);
    
elseif status==0
    disp('The RESULT DISPLAY module did not complete successfully. Aborting.');
    return
else
    disp('The RESULT DISPLAY module was not run.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   SAVE fsmParam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mark modules up-to-date
fsmParam=fsmMarkModuleUpToDate(fsmParam);

% Add date and time of the last run
dateFinished=datestr(now,0);
fsmParam.specific.lastRun=dateFinished; 

% Save fsmParam to disk
eval(['save ',fsmParam.main.path,filesep,'fsmParam.mat fsmParam']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   STOP TIMER
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Store elapsed time
elapsedTime=toc;

% Write to console
fprintf(1,'All done! Elapsed time: %d seconds.\n',round(elapsedTime));

% Write as well to 'parameters.txt'
fid=fopen('log.txt','a');

% Could the file be opened successfully?
if fid==-1
	error('Couldn''t open the file.');
end

% Write to file
fprintf(fid,'\n\n****************************************************************************************\n');
fprintf(fid,'*                                                                                      *\n');
fprintf(fid,'*     Calculation finished: %s                                       *\n',dateFinished);
fprintf(fid,'*                                                                                      *\n');
fprintf(fid,'****************************************************************************************\n\n\n');

% Close the file
fclose(fid);

% Save image information for other packages
fsmImages.firstIndex = fsmParam.specific.firstIndex;
fsmImages.firstName  = fsmParam.specific.fileList(1,:);
fsmImages.lastIndex  = fsmParam.specific.lastIndex;
fsmImages.lastName   = fsmParam.specific.fileList(end,:);
save fsmImages.mat fsmImages;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   CHANGE BACK TO PREVIOUS DIRECTORY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(olddir);


