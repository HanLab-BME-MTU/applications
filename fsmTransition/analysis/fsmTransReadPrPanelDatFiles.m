function [normals,profiles,protrusions,edgePixels,segments]=fsmTransReadPrPanelDatFiles(edgePixelsFileName,normalsFileName,toggleProfiles,protFileName,segmFileName,vLength)
% fsmTransReadPrPanelDatFiles imports prPanel .dat files into MATLAB matrices/structures
%
% SYNOPSIS   [normals,profiles,protrusions,edgePixels,segments]=fsmTransReadPrPanelDatFiles(imgFileName,edgePixelsFileName,normalsFileName,protFileName,vLength)
%
% INPUT      edgePixelsFileName: string containing the file name (with full path) of the pixel_edge.dat file
%                                Pass edgePixelsFileName=[] to select the image through a dialog.
%                                Pass edgePixelsFileName=-1 to skip it.
%            normalsFileName   : string containing the file name (with full path) of the av_normal.dat file
%                                Pass normalsFileName=[] to select the image through a dialog.
%                                Pass normalsFileName=-1 to skip it.
%            toggleProfiles    : [ 0 | 1 ] Turns off | on the creation of profiles (see OUTPUT).
%            protFileName      : string containing the file name (with full path) of the av_prot.dat file
%                                Pass protFileName=[] to select the image through a dialog.
%                                Pass protFileName=-1 to skip it.
%            vLength           : length of the profile vectors. Optional, default = 100;
%
% OUTPUT     normals           : (n x 4 x k) matrix, where n is the number of segments along the edge (as set in prPanel) and
%                                k is the number of processed images. [y0 x0 y x]n are the coordinates of the normal vectors to
%                                the extracted edge.
%            profiles          : (n x 4 x k) matrix; n segments, k images. Vectors perpendicular to the edge, pointing within the 
%                                object (opposite to the normals) of length vLength
%            protrusions       : (n x 4 x k) matrix; n segments, k images. Protrusion vectors calculated at the n segments.
%            edgePixels        : structure(1:k).x : x coordinates of pixel edges 
%                                              .y : y coordinates of pixel edges
%            segments          : coordinates of segments along the edge
%
% Aaron Ponti, 04/26/2004

% Check input parameters
if nargin<6 | nargin>7
    error('6 or 7 input parameters exepected');
end

% Set default value for vLength if needed
if nargin==6
    vLength=100;
end

% Toggle flags
% Pixel edge
if isempty(edgePixelsFileName) 
    PIXELEDGE=1; % The user wants to pick the file
else
    if edgePixelsFileName==-1
        PIXELEDGE=0; % The user does not want this file to be imported (-1)
    else
        PIXELEDGE=1;
        if exist(edgePixelsFileName)~=2
            edgePixelsFileName=[]; % The user's specified does not exist. He will have to pick the file
        end
    end
end
% Normals
if isempty(normalsFileName) 
    NORMALS=1; % The user wants to pick the file
else
    if normalsFileName==-1
        NORMALS=0; % The user does not want this file to be imported (-1)
    else
        NORMALS=1;
        if exist(normalsFileName)~=2
            normalsFileName=[]; % The user's specified does not exist. He will have to pick the file
        end
    end
end
% Protrusion
if isempty(protFileName) 
    PROTRUSION=1; % The user wants to pick the file
else
    if protFileName==-1
        PROTRUSION=0; % The user does not want this file to be imported (-1)
    else
        PROTRUSION=1;
        if exist(protFileName)~=2
            protFileName=[]; % The user's specified does not exist. He will have to pick the file
            
        end
    end
end
% Profiles
if NORMALS==0
    PROFILES=0; % If normals are not loaded, profiles cannot be created
else
    PROFILES=toggleProfiles;
end
if isempty(segmFileName) 
    SEGMENTS=1; % The user wants to pick the file
else
    if segmFileName==-1
        SEGMENTS=0; % The user does not want this file to be imported (-1)
    else
        SEGMENTS=1;
        if exist(segmFileName)~=2
            segmFileName=[]; % The user's specified does not exist. He will have to pick the file
        end
    end
end


% Initialize outputs
normals=[]; profiles=[]; protrusions=[]; edgePixels=[]; segments=[];

% Check that at least one field is active
if sum([PIXELEDGE NORMALS PROTRUSION SEGMENTS])==0
    disp('Nothing to do here.');
    return
end

% Initialize wiatbar
hWait=waitbar(0,'Please wait...');
nSteps=4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CELL EDGE COORDINATES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PIXELEDGE~=0
       
    % Do we need to ask the user to pick a file?
    if isempty(edgePixelsFileName)
        
        % Load pixel_edge.dat
        [fName,dirName] = uigetfile(...
            {'pixel_edge.dat;','pixel_edge.dat';
            '*.*','All Files (*.*)'},...
            'Select pixel_edge.dat');
        if ~(isa(fName,'char') & isa(dirName,'char'))
            return 
        end
    else
        
        [dirName,body,no,ext]=getFilenameBody(edgePixelsFileName);
        fName=[body,no,ext];
    end
    
    % Open pixel_edge.dat file
    [fid0,message]=fopen([dirName,filesep,fName],'r');
    if fid0==-1
        error('Could not open file');
    end

    % Create structure edgePixels
    edgePixels=struct('y',[],'x',[]);
    
    linesPerFrame=3;
    frame=1;
    lineNo=1;
    while not(feof(fid0)) 
        
        % Check where we are
        current=mod(lineNo,linesPerFrame);
        if current==0
            current=linesPerFrame;
        end
        
        % Read new line from both files
        tline=fgetl(fid0);
        
        % Break if end of file
        if ~ischar(tline)
            break
        end
        
        if isempty(tline)
            % Jump to next line
            continue
        end
        
        % Write the coordinates
        switch current
            case 2, 
                edgePixels(frame).x=str2num(tline)'; 
            case 3, 
                edgePixels(frame).y=str2num(tline)';
            otherwise,
        end
        
        % Next line
        lineNo=lineNo+1;
        
        % Update frame
        if current==linesPerFrame
            frame=frame+1;
        end
        
    end
    fclose(fid0);
   
end
% Update waitbar
waitbar(1/nSteps,hWait);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NORMALS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if NORMALS~=0
    
    % Do we need to ask the user to pick a file?
    if isempty(normalsFileName)
        
        % Load av-normal.dat
        [fName,dirName] = uigetfile(...
            {'av_normal.dat;','av_normal.dat';
            '*.*','All Files (*.*)'},...
            'Select av_normal.dat');
        if ~(isa(fName,'char') & isa(dirName,'char'))
            return 
        end
    else
        
        [dirName,body,no,ext]=getFilenameBody(normalsFileName);
        fName=[body,no,ext];
    end
    
    % Open av-normal.dat file
    [fid,message]=fopen([dirName,filesep,fName],'r');
    if fid==-1
        error('Could not open file');
    end

    linesPerFrame=6;
    frame=1;
    lineNo=1;
    while not(feof(fid)) 
        
        % Check where we are
        current=mod(lineNo,linesPerFrame);
        if current==0
            current=linesPerFrame;
        end
        
        % Read new line from both files
        tline=fgetl(fid);
        
        % Break if end of file
        if ~ischar(tline)
            break
        end
        
        if isempty(tline)
            % Jump to next line
            continue
        end
        
        % Write the coordinates
        switch current
            case 3, 
                normals(:,2,frame)=str2num(tline)'; 
            case 4, 
                normals(:,1,frame)=str2num(tline)';
            case 5, 
                normals(:,4,frame)=str2num(tline)'+normals(:,2,frame);
            case 6, 
                normals(:,3,frame)=str2num(tline)'+normals(:,1,frame);
            otherwise,
        end
        
        % Next line
        lineNo=lineNo+1;
        
        % Update frame
        if current==linesPerFrame
            frame=frame+1;
        end
        
    end
    fclose(fid);
   
    if PROFILES==1
        % Calculate profiles
        normalVectors=[normals(:,3,:)-normals(:,1,:) normals(:,4,:)-normals(:,2,:)];
        normalLengths=sqrt(normalVectors(:,1,:).^2+normalVectors(:,2,:).^2);
        profileLengths(:,1,:)=vLength.*normalVectors(:,1,:)./normalLengths;
        profileLengths(:,2,:)=vLength.*normalVectors(:,2,:)./normalLengths;
        profiles=[normals(:,1,:) normals(:,2,:) normals(:,1,:)-profileLengths(:,1,:) normals(:,2,:)-profileLengths(:,2,:)];
    end
    
end

% Update waitbar
waitbar(2/nSteps,hWait);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROTRUSION VECTORS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PROTRUSION~=0

    % Do we need to ask the user to pick a file?
    if isempty(protFileName)
        
        % Load av-prot.dat
        [fName,dirName] = uigetfile(...
            {'av_prot.dat;','av_prot.dat';
            '*.*','All Files (*.*)'},...
            'Select av_prot.dat');
        if ~(isa(fName,'char') & isa(dirName,'char'))
            return 
        end
    else
        
        [dirName,body,no,ext]=getFilenameBody(protFileName);
        fName=[body,no,ext];
        
    end

    
    % Open av-prot.dat file
    [fid2,message]=fopen([dirName,filesep,fName],'r');
    if fid2==-1
        error('Could not open file');
    end
    
    linesPerFrame=6;
    frame=1;
    lineNo=1;
    while not(feof(fid2)) 
        
        % Check where we are
        current=mod(lineNo,linesPerFrame);
        if current==0
            current=linesPerFrame;
        end
        
        % Read new line from both files
        tline2=fgetl(fid2);
        
        % Break if end of file
        if ~ischar(tline2)
            break
        end
        
        if isempty(tline2)
            % Jump to next line
            continue
        end
        
        % Write the coordinates
        switch current
            case 3, 
                protrusions(:,2,frame)=str2num(tline2)';
            case 4, 
                protrusions(:,1,frame)=str2num(tline2)';
            case 5, 
                protrusions(:,4,frame)=str2num(tline2)'+protrusions(:,2,frame);
            case 6, 
                protrusions(:,3,frame)=str2num(tline2)'+protrusions(:,1,frame);
            otherwise,
        end
        
        % Next line
        lineNo=lineNo+1;
        
        % Update frame
        if current==linesPerFrame
            frame=frame+1;
        end
        
    end
    fclose(fid2);
    
end

% Update waitbar
waitbar(3/nSteps,hWait);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SEGMENTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if SEGMENTS~=0
    
    % Load regions
    if isempty(segmFileName)
        
        % Load s_mask###.dat
        [fName,dirName] = uigetfile(...
            {'*.dat;','*.dat';
            '*.*','All Files (*.*)'},...
            'Select s_mask###.dat');
        if ~(isa(fName,'char') & isa(dirName,'char'))
            return 
        end
    else
        
        [dirName,body,no,ext]=getFilenameBody(segmFileName);
        fName=[body,no,ext];
    end
    
    outFileList2=getFileStackNames([dirName,filesep,fName]);
    
    segments=struct('y',[],'x',[]);
    
    for i=1:length(outFileList2)
        
        % Open s_mask###.dat file
        [fid3,message]=fopen(char(outFileList2(i)),'r');
        if fid3==-1
            error('Could not open file');
        end
        
        linesPerFrame=3;
        frame=1;
        lineNo=1;
        while not(feof(fid3))    
            
            % Check where we are
            current=mod(lineNo,linesPerFrame);
            if current==0
                current=linesPerFrame;
            end
            
            % Read new line from both files
            tline3=fgetl(fid3);
            
            % Break if end of file
            if ~ischar(tline3)
                break
            end
            
            if isempty(tline3)
                % Jump to next line
                continue
            end
            
            % Write the coordinates
            switch current
                case 2,
                    polyX=str2num(tline3)';
                    segments(i).x(frame,1:length(polyX))=polyX;
                case 3, 
                    polyY=str2num(tline3)';
                    segments(i).y(frame,1:length(polyY))=polyY;            
                otherwise,
            end
            
            % Next line
            lineNo=lineNo+1;
            
            % Update frame
            if current==linesPerFrame
                frame=frame+1;
            end
            
            %     pause;
            
        end
        fclose(fid3);
        
    end
    
end

% Update waitbar
waitbar(4/nSteps,hWait);

% Close waitbar
close(hWait);