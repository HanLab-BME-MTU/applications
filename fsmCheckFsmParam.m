function [fsmParam,status]=fsmCheckFsmParam(fsmParam)
% fsmCheckFsmParam checks that the 'specific' fields of fsmParam.mat are valid  
%
% SYNOPSIS      [fsmParam,status]=fsmCheckFsmParam(fsmParam)
%
% INPUT         fsmParam :   general parameter structure for fsm
%                            Type help fsmGetParamDflts for more info on fsmParam
%
% OUTPUT        fsmParam :   (corrected) fsmParam 
%               status   :   success flag. 0 : the function did not complete successfully
%                                          1 : the function completed successfully
%                            
% DEPENDENCES              
%
% Aaron Ponti, April 26th, 2004


if nargin~=1
    error('One input parameter expected');
end

% Default status is 0 (unsuccessful)
status=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% IF fsmParam.mat IS NOT THE DEFAULT FILE AND CONTAIN ALREADY SOME SPECIFIC
%   EXPERIMENT SETTINGS, CHECK THAT THEY ARE VALID - ACTUALLLY, ONLY THE LIST
%   OF IMAGES IS CHECKED FOR VALIDITY HERE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(fsmParam.specific.fileList) % Already analyzed
    if ~exist(fsmParam.specific.fileList(1,:))
        infoMsg=['The stored list of images is not valid. The first image [',fsmParam.specific.fileList(1,:),'] does not exist. You will now be asked to locate this image.'];
        uiwait(msgbox(infoMsg,'Info','modal'));
        
        % Create title for open dialog
        [path,body,noOriginal,ext]=getFilenameBody(fsmParam.specific.fileList(1,:));
        
        title=['Please locate ', body,noOriginal,ext];
        
        valid=0;
        while valid==0
            
            % The user must select the first image of the stack 
            [fName,dirName] = uigetfile(...
                {'*.tif;*.tiff;*.jpg;*.jpeg','Image Files (*.tif,*.tiff,*.jpg,*.jpeg)';
                '*.tif','TIF files (*.tif)'
                '*.tiff','TIFF files (*.tiff)'
                '*.jpg;','JPG files (*.jpg)'
                '*.jpeg;','JPEG files (*.jpeg)'
                '*.*','All Files (*.*)'},...
                title);
            if(isa(fName,'char') & isa(dirName,'char'))
                
                if exist([dirName,fName])==2
                    
                    % Check that the frame number matches
                    [path,body,noNew,ext]=getFilenameBody([dirName,fName]);
                    if strcmp(noOriginal,noNew)==0
                        
                        % Prepare error message
                        errorMsg='The selected frame number does not match.';
                        
                    else
                        
                        % Recover all file names from the stack
                        outFileList=getFileStackNames([dirName,fName]);
                        
                        % Take only the number of images specified in fsmParam.specific.imageNumber
                        if fsmParam.specific.imageNumber<=length(outFileList)
                            
                            outFileList=outFileList(1:fsmParam.specific.imageNumber);
                            
                            % Convert filelist to a matrix of strings
                            outFileList=char(outFileList);
                            
                            % Replace the current list of files in fsmParam.specific.fileList
                            fsmParam.specific.fileList=outFileList;
                            
                            % Set flags to success
                            valid=1; status=1;
                            
                        else
                            
                            % Prepare error message
                            errorMsg=['Last time, ',num2str(fsmParam.specific.imageNumber),' were processed. This stack contains only ',num2str(length(outFileList)),' images.'];
                            
                        end
                        
                    end
                else
                    
                    % Prepare error message
                    errorMsg='The specified file does not exist.'; % This should not happen
                    
                end
                
            else

                % Prepare error message
                errorMsg='No file specified.';

            end
            
            if valid==0 % Still not fixed
                
                choice=questdlg([errorMsg,' Try again?'], ...
                    'Error', ...
                    'Yes','No','Yes');
                switch choice
                    case 'Yes', valid=0;
                    case 'No',  valid=1; status=0; % This informs fsmMain to quit
                    otherwise
                        error('Wrong selection');
                end
                
            end
            
        end
            
    else
        
        % The images exist, inform fsmMain that it can continue safely.
        status=1;
        
    end
    
end
