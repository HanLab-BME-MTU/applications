function fsmWriteParamsToTextFile(fsmParam,fsmExpParam,userPath)
% fsmWriteParamsToTextFile writes the parameters from fsmParam into a text file
%
% SYNOPSIS      fsmWriteParamsToTextFile(fsmParam,fsmExpParam,userPath)
%
% INPUT         fsmParam    : parameter structure
%                             (type 'help fsmGetparamDflts' for more detail)
%               fsmExpParam : structure of noise parameters as created by fsmGuiScanDataBase
%               userPath    : work path selected in the user interface
%
% OUTPUT        None - a text file called 'parameters.txt' is written in 'userPath'.
%
% DEPENDENCES   fsmMain uses { }
%               fsmMain is used by { fsmGuiMain }
%
% Aaron Ponti, May 5th, 2003

% Check input parameters
if nargin~=3
    error('Three parameters expected');
end

% Check the analysis has been run, otherwise return
if isempty(fsmParam.specific.fileList)

    % The software has been interrupted - nothing to write
    return;
    
end

% Change to the current path
oldPath=cd;

cd(userPath);

% Open/create text file 
fid=fopen('parameters.txt','a+');
if fid==-1
    error('Cannot open/create file.');
end

% Write date
fprintf(fid,'Created: %s\n',datestr(now));

% Write title
fprintf(fid,'\nParameters\n----------\n\n');

% Start writing fields

% MAIN & SPECIFIC
fprintf(fid,'[ EXPERIMENT INFO ]\n\n');
fprintf(fid,'Work path                 : %s\n',fsmParam.main.path);
fprintf(fid,'Number of images          : %d\n',fsmParam.specific.imageNumber); % This must be read from .specific
fprintf(fid,'Image size                : %d x %d\n',fsmParam.specific.imgSize(1),fsmParam.specific.imgSize(2)); % This must be read from .specific
fprintf(fid,'First image name and path : %s\n',fsmParam.specific.fileList(1,:)); % This must be read from .specific
fprintf(fid,'Bit depth                 : %d\n',log2(fsmParam.main.normMax+1));
qtile=[1.15 1.29 1.45 1.645 1.96 2.58];
pcent=[75 80 85 90 95 99];
quantile=fsmParam.main.noiseParam(5);
percent=pcent(find(qtile==quantile));
if isempty(percent)
    percent='user-defined confidence level';
else
    percent=[num2str(percent),'%'];
end
fprintf(fid,'Quantile                  : %1.2f (%s)\n',quantile,percent);
fprintf(fid,'Noise parameters          : [ sDN = %1.8f; beta = %1.8f; I0 = %1.8f ]\n',fsmParam.main.noiseParam(2),fsmParam.main.noiseParam(3),fsmParam.main.noiseParam(4));

% PREP
fprintf(fid,'Experiment name           : %s\n',fsmExpParam(fsmParam.main.noiseParam(7)-1).label); % This belongs to MAIN
fprintf(fid,'Experiment description    : %s\n',fsmExpParam(fsmParam.main.noiseParam(7)-1).description);

fprintf(fid,'\n[ MODULES ]\n');

fprintf(fid,'\n> PREP \n');
if fsmParam.prep.enable==1
    
    fprintf(fid,'Module ''prep''             : run\n');
    switch fsmParam.prep.pstSpeckles
        case 1, fprintf(fid,'Speckle order             : Primary only\n');
        case 2, fprintf(fid,'Speckle order             : %d-ary (Minimum increase: %1.3f)\n',fsmParam.prep.paramSpeckles(1),fsmParam.prep.paramSpeckles(2));
        case 3, fprintf(fid,'Speckle order             : Scale space (sigma: &1.3f)\n',fsmParam.prep.paramSpeckles(1));
        otherwise
            error('Wrong value for fsmParam.prep.pstSpeckles');
    end
    switch fsmParam.prep.enhTriang
        case 1, fprintf(fid,'Enhanced triangulation    : Yes\n');
        case 0, fprintf(fid,'Enhanced triangulation    : No\n');
        otherwise
            error('Wrong value for fsmParam.prep.enhTriang');
    end
    switch fsmParam.prep.autoPolygon
        case 1, fprintf(fid,'Automatic edge detection  : Yes\n');
        case 0, fprintf(fid,'Automatic edge detection  : No\n');
        otherwise
            error('Wrong value for fsmParam.prep.autoPolygon');
    end
    fprintf(fid,'Gauss ratio               : %1.2f\n',fsmParam.prep.gaussRatio);
    
elseif fsmParam.prep.enable==0
    
    fprintf(fid,'Module ''prep''             : not run\n');
    
else

    error('Wrong value for smParam.prep.enable');
    
end

% TRACK
fprintf(fid,'\n> TRACK \n');

if fsmParam.track.enable==1
    
    fprintf(fid,'Module ''track''            : run\n');
    switch fsmParam.track.tracker
        case 1, fprintf(fid,'Selected tracker          : Brownian Motion Tracker + Neural Network\n');
        case 2, fprintf(fid,'Selected tracker          : Enhanced Brownian Motion Tracker\n');
        case 3, fprintf(fid,'Selected tracker          : Flow tracker (3 frames)\n');
        otherwise
            error('Wrong value for fsmParam.prep.tracker');
    end
    switch fsmParam.track.enhanced
        case 1, 
            fprintf(fid,'Hierarchical tracking     : Yes\n');    
            switch fsmParam.track.grid
                case 1, fprintf(fid,'-> Grid enabled           : Yes\n');
                case 0, fprintf(fid,'-> Grid enabled           : No\n');
                otherwise
                    error('Wrong value for fsmParam.prep.grid');
            end
        case 0, fprintf(fid,'Hierarchical tracking     : No\n');
        otherwise
            error('Wrong value for fsmParam.prep.enhanced');
    end

    fprintf(fid,'Search radius (pixels)    : %d\n',fsmParam.track.threshold);
    
elseif fsmParam.track.enable==0
    
    fprintf(fid,'Module ''track''            : not run\n');
    
else
    error('Wrong value for fsmParam.track.enable');
    
end

% BUILD
fprintf(fid,'\n> BUILD \n');
switch fsmParam.build.enable
    case 1, fprintf(fid,'Module ''build''            : run\n');
    case 0, fprintf(fid,'Module ''build''            : not run\n');
    otherwise
        error('Wrong value for fsmParam.build.enable');
end

% KIN
fprintf(fid,'\n> KIN \n');
if fsmParam.kin.enable==1
    
    fprintf(fid,'Module ''kin''              : run\n');
    switch fsmParam.kin.bleachRed
        case 0,        fprintf(fid,'Bleaching reductiom       : Off\n');;
        case 7.25e-5,  fprintf(fid,'Bleaching reductiom       : 1x\n');
        case 1.45e-4,  fprintf(fid,'Bleaching reductiom       : 2x\n');
        case 2.175e-4, fprintf(fid,'Bleaching reductiom       : 3x\n');
        otherwise 
            error('Wrong value for fsmParam.kin.bleachRed');
    end
    
elseif fsmParam.kin.enable==0
    
    fprintf(fid,'Module ''kin''              : not run\n');
    
else
    
    error('Wrong value for fsmParam.kin.enable');
    
end


% DISP
fprintf(fid,'\n> DISP \n');
switch fsmParam.disp.enable
    case 1, fprintf(fid,'Module ''disp''             : run\n');
    case 0, fprintf(fid,'Module ''disp''             : not run\n');
    otherwise 
        error('Wrong value for fsmParam.disp.enable');
end

% Add separation
fprintf(fid,'\n##################################################################\n');

% Close file
fclose(fid);

% Change back to original path
cd(oldPath);
