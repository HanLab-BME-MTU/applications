function fsmExpParam=fsmGuiScanDataBase(fileName)
% fsmGuiScanDataBase reads the experimental parameters from a file and fills the structure fsmExpParam
%
% SYNOPSIS      fsmExpParam=fsmGuiScanDataBase(fileName)
%
% INPUT         filename    :   file name with path of the parameter file
%
% OUTPUT        fsmExpParam :   structure array containing the experimental parameters.
%                               It contains following fields:
%                                           (fsmExpParam).label
%                                                        .bitdepth
%                                                        .noiseparams
%                                                        .gaussratio
% DEPENDENCES   
%               
% Aaron Ponti, March 27th, 2003

% Initialize structure fsmExpParam to store experiment parameters
fsmExpParam=struct('label','',...
    'description','',...
    'bitDepth',0,...
    'noiseParams',[0 0 0],...
    'gaussRatio',0);

% Experiment counter
counter=0;

% Set parameter flags to zero
L_=0; D_=0; B_=0; N_=0; G_=0;

% Open file
fid=fopen(fileName);

% Check whether the file exists
if fid==-1
    error('Could not find the file.');
end

% Scan the file and fill fsmExpParam
while not(feof(fid))
    
    % Read new line
    tline=fgetl(fid);
    
    % Break if end of file
    if ~ischar(tline)
        break
    end
    
    if isempty(tline)
        % Jump to next line
        continue
    end
    
    if tline(1)=='%'
        % This is a comment - ignore it
        continue
    end
        
    % This is the beginning of an experiment (the sign '#')
    if tline(1)=='#'
        
        if counter==0 
            
            % This is the first experiment found
            counter=counter+1;
            
            % Turn off bits
            L_=0; D_=0; B_=0; N_=0; G_=0;
            
        else

            % Save parameters for previous experiment
            if L_==1 & D_==1 & B_==1 & N_==1 & G_==1
                
                % All parameters have been read for this experiment
                fsmExpParam(counter).label=label;
                fsmExpParam(counter).description=description;
                fsmExpParam(counter).bitDepth=bitdepth;
                fsmExpParam(counter).noiseParams=noiseparams;
                fsmExpParam(counter).quantile=vQuantile;
                fsmExpParam(counter).gaussRatio=gaussratio;
              
            end
            
            % Start collecting data for the next experiment
            counter=counter+1;
            
            % Turn off bits
            L_=0; D_=0; B_=0; N_=0; G_=0;
            
            % Clear variables
            label=[]; bitdepth=[]; noiseparams=[]; gaussratio=[]; vQuantile=[];

        end

        % Store experiments
        
    elseif findstr(tline,'LABEL')
        
        if L_==1
            % More than one line LABEL for this experiment
            error('The file format does not follow the specifications.');
        end
        
        % Read label
        label=readString(tline,fid);
        
        % Set bit for label to 1
        L_=1;
        
    elseif findstr(tline,'DESCRIPTION')
        
        if D_==1
            % More than one line DESCRIPTION for this experiment
            error('The file format does not follow the specifications.');
        end
        
        % Read bit depth
        description=readString(tline,fid);
        
        % Set bit for bitdepth to 1
        D_=1;
        
    elseif findstr(tline,'BIT DEPTH')
        
        if B_==1
            % More than one line BIT DEPTH for this experiment
            error('The file format does not follow the specifications.');
        end
        
        % Read bit depth
        bitdepth=str2num(readString(tline,fid));
        
        % Set bit for bitdepth to 1
        B_=1;

    elseif findstr(tline,'NOISE PARAMS')

        if N_==1
            % More than one line NOISE PARAMS for this experiment
            error('The file format does not follow the specifications.');
        end

        % Read noise parameters
        noiseparams=str2num(readString(tline,fid));
        if length(noiseparams)==3
            vQuantile=0;
        elseif length(noiseparams)==4
            vQuantile=noiseparams(4);
            noiseparams=noiseparams(1:3);            
        end
        
        % Set bit for bitdepth to 1
        N_=1;

    elseif findstr(tline,'GAUSS RATIO')

        if G_==1
            % More than one line GAUSS RATIO for this experiment
            error('The file format does not follow the specifications.');
        end

        % Read noise parameters
        gaussratio=str2num(readString(tline,fid));
        
        % Set bit for bitdepth to 1
        G_=1;

    else
        
        % A line beginning with any other char - ignore the line
        
    end

end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function string=readString(tline,fid)

% Find the string
indx=find(tline=='"');
if length(indx)~=2
    fclose(fid);
    error('The file format does not follow the specifications.');
end
        
% Read string
string=tline(indx(1)+1:indx(2)-1);
