function speckleArray=fsmBuildSaveSpeckleArray(cM,outFileList,xmin,xmax,firstIndex,strg,noiseParams,threshold,factors,W,userPath,sigma)
% fsmBuildSaveSpeckleArray fills the speckleArray structure containing all speckle information extracted from the analysis of an FSM stack
%
% SYNOPSIS   speckleArray=fsmBuildSaveSpeckleArray(cM,outFileList,xmin,xmax,firstIndex,strg,noiseParams,threshold,factors,userPath,sigma)
%
% INPUT      cM            : magic position matrix (mpm) as returned by the tracker module
%            outFileList   : string matrix containing all image file names.
%                            It has the form:
%                                                   |'img/img01.tif'|
%                                                   |'img/img02.tif'|
%                                     outFileList = |'img/img03.tif'|
%                                                   |'img/img04.tif'|
%                                                   |'img/img05.tif'|                                
%
%            xmin          : lower boundary for image normalization
%            xmax          : upper boundary for image normalization
%            firstIndex    : index of the corresponding first image
%            strg          : format string to correctly format the numeric suffix of saved files
%            noiseParams   : noise parameters for speckle significance test
%            threshold     : radius of the region considered when looking for matching speckles
%            factors       : vector [1x(image number)] containing the intensity correction
%                            factors (between 0 and 1) - this is TURNED OFF in FSM because may bias
%                            the result - pass a vector full of 1's to turn it off.
%            W             : [y0 x0 y x], sets the upper-left (y0,x0) and lower-right (y,x) 
%                            corners if the analyzed region is smaller than the image 
%                            (e.g. when a grid is created with framework).
%                            Otherwise, set W=[1 1;size(image,1) size(image,2)]
%            userPath      : work path
%            sigma         : sigma for image low-pass filtering (optional, default = 1)
%
% OUTPUT     speckleArray  : structure containing all speckle information extracted from an image stack
%
% fsmBuildSaveSpeckleArray is used by { fsmBuildMain }
% fsmBuildSaveSpeckleArray uses { }

if nargin<11 || nargin>12
    error('Wrong number of input arguments');
end

% Sigma has a default value of 1
if nargin==11
    sigma=1;
end

tp=size(cM,2)/2; % Number of time points
tm=size(cM,3);   % Number of subwindows
sp=size(cM,1);   % Number of speckle lines

% Initialize waitbar 
h=waitbar(0,'Indexing speckles...');
tot=sp*tp;

% Calculate some longer step size (not to spend too much time updating the waitbar)
step=0.5*10^fix(log10(tot)-1); if step<1, step=1; end

%
% Index speckles and non-speckles
%

x=1:2:2*tp-1;

% Allocate memory for 'table'
table=zeros(sp*tp,4);

pos4=0;                               % This marks a speckle still existing in the last time point
cMPos=0;                              % Position in cM
sAPos=0;                              % Position in speckleArray
for c1=1:sp
    for c2=1:tp
        cMPos=cMPos+1;                  % Keep info of the position in cM  
        currentSpeckle=cM(c1,x(c2));    % Read x coordinate of current speckle
        
        %
        % We are at the beginning of the movie
        %
        
        if c2==1
            
            if currentSpeckle~=0         % Here we have a speckle
                pos4=1;                  % We have to mark this because we have a speckle already existing at the first time point                               
                if cM(c1,x(c2+1))~=0      % Also next one is speckle
                    sAPos=sAPos+1;         % Increase also the position in the struct by 1
                    pos1=0;         
                    pos2=sAPos;
                    pos3=0;
                else                      % Next speckle is zero: here we have to introduce a post-death non-speckle
                    sAPos=sAPos+1;         % This is the position of current speckle
                    pos1=0;
                    pos2=sAPos;
                    sAPos=sAPos+1;         % This will be the position of the "post-death" non-speckle
                    pos3=sAPos;
                end
            else                         % We have a non-speckle
                pos1=-1;
                pos2=-1;
                pos3=-1;
            end
            
            %
            % We are at the end of the movie
            %
            
        elseif c2==tp
            
            if currentSpeckle~=0         % Here we have a speckle
                pos4=2;                  % We have to cut here because we are at the last time point with a still existing speckle
                if cM(c1,x(c2-1))==0      % Last one was a non-speckle
                    sAPos=sAPos+1;
                    pos1=sAPos;            % This is the "pre-birth" non-speckle
                    sAPos=sAPos+1;
                    pos2=sAPos;            % This is the speckle
                    pos3=0; 
                else                      % Last one was also a speckle
                    sAPos=sAPos+1;         % Increase also the position in the struct by 1
                    pos1=0;          
                    pos2=sAPos;
                    pos3=0;
                end
            else                         % If we are at the last time-point (and we have a speckle)
                pos1=-1;
                pos2=-1;  
                pos3=-1;
            end
            
            %
            % We are in the central time-points
            %   
            
        else
            
            if currentSpeckle~=0          % Here we have a speckle
                
                if cM(c1,x(c2-1))~=0 && cM(c1,x(c2+1))~=0   % Last and next ones are speckle
                    pos1=0;
                    sAPos=sAPos+1;          % This is the speckle
                    pos2=sAPos;
                    pos3=0;
                end
                
                if cM(c1,x(c2-1))==0 && cM(c1,x(c2+1))~=0   % Last one was not a speckle; the next is
                    sAPos=sAPos+1;
                    pos1=sAPos;             % This is the "pre-birth" non-speckle
                    sAPos=sAPos+1;
                    pos2=sAPos;             % This is the speckle
                    pos3=0;
                end
                
                if cM(c1,x(c2-1))~=0 && cM(c1,x(c2+1))==0   % Last one was a speckle; the next is not
                    pos1=0;
                    sAPos=sAPos+1;          % This is the speckle
                    pos2=sAPos;
                    sAPos=sAPos+1;          % This is the "post-death" non-speckle
                    pos3=sAPos;
                end
                
                if cM(c1,x(c2-1))==0 && cM(c1,x(c2+1))==0   % Both are not speckles
                    sAPos=sAPos+1;
                    pos1=sAPos;             % This is the "pre-birth" non-speckle 
                    sAPos=sAPos+1;
                    pos2=sAPos;             % This is the speckle
                    sAPos=sAPos+1;
                    pos3=sAPos;             % This is the "post-death" non-speckle
                end
                
            else                         % We have a non-speckle
                pos1=-1;
                pos2=-1;
                pos3=-1;
            end
        end
        % Write correspondence
        table(cMPos,1:4)=[pos1 pos2 pos3 pos4];
        % Reset pos4;
        pos4=0;
        
        % Update waitbar if needed
        if step>1
            if mod((c1*c2),step)==1
                waitbar((c1*c2)/tot,h);
            end
        else
            waitbar((c1*c2)/tot,h);
        end        
        
    end
end

close(h);

%
% Write speckle info into struct
%

% Allocate memory for speckleArray
lengthSpeckleArray=max(table(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize empty speckle array
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

speckleArray=struct(...
    'timepoint',   uint16(zeros(lengthSpeckleArray,1)),...  % Unsigned integer 16: 1 - 65535 (2 bytes)
    'spPos',       uint16(zeros(lengthSpeckleArray,2)),...  % Unsigned integer 16 
    'bgPos1',      uint16(zeros(lengthSpeckleArray,2)),...  % Unsigned integer 16
    'bgPos2',      uint16(zeros(lengthSpeckleArray,2)),...  % Unsigned integer 16
    'bgPos3',      uint16(zeros(lengthSpeckleArray,2)),...  % Unsigned integer 16
    'intensity',   zeros(lengthSpeckleArray,1),...          % Double (8 bytes)
    'background',  zeros(lengthSpeckleArray,1),...          % Double
    'deltaI',      zeros(lengthSpeckleArray,1),...          % Double
    'deltaICrit',  zeros(lengthSpeckleArray,1),...          % Double
    'sigmaSp',     zeros(lengthSpeckleArray,1),...          % Double 
    'sigmaBg',     zeros(lengthSpeckleArray,1),...          % Double
    'status',      char(zeros(lengthSpeckleArray,1)),...    % Char (1 byte)
    'speckleType', uint8(zeros(lengthSpeckleArray,1)),...   % Unsigned integer 8: 1 - 255 (1 byte)
    'score',       zeros(lengthSpeckleArray,1),...          % Double 
    'activity',    int8(zeros(lengthSpeckleArray,1)),...    % Integer 8: -128 - 127 (1 byte)
    'lmEvent',     false(lengthSpeckleArray,1));            % Logical: 0 | 1 (1 byte)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize waitbar 
h=waitbar(0,'Populating speckleArray structure...');
tot=size(cM,1)*0.5*size(cM,2)*size(cM,3);

% Calculate some longer step size (not to spend too much time updating the waitbar)
step=0.5*10^fix(log10(tot)-2); if step<1, step=1; end

for c1=1:tp                                    % Cycle through timepoints
    
    % Current index
    currentIndex=c1+firstIndex-1;
    
    % Index string
    indxStr=sprintf(strg,currentIndex);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % gapList
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Load gapList
    eval(['load ',userPath,filesep,'gapList',filesep,'gapList',indxStr,'.mat;']); % gapList
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % cands
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Current
    eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'.mat;']); % Cands workspace
    candsS=cands; % Structure as is
    Lmax=[cands.Lmax]; 
    allCandsPos=reshape(Lmax,2,length(Lmax)/2)'; % All local maximum positions
    if ~isempty(cands)
        insCands=cands(find([cands.status]==0)); % Cands structure containing only insignificant loc max
        clear Lmax;
        Lmax=[insCands.Lmax]; 
        insCandsPos=reshape(Lmax,2,length(Lmax)/2)'; % Insignificant loc max positions
    else
        insCands=[];
        insCandsPos=[];
    end
    clear cands;
    
    % Previous
    if c1>1
        indxStr=sprintf(strg,currentIndex-1);
        eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'.mat;']); % Cands workspace
        candsB=cands; % Structure as is
        Lmax=[cands.Lmax]; 
        allCandsBPos=reshape(Lmax,2,length(Lmax)/2)'; % All local maximum positions
        if ~isempty(cands)
            insCandsB=cands(find([cands.status]==0)); % Cands structure containing only insignificant loc max
            clear Lmax;
            Lmax=[insCandsB.Lmax]; 
            insCandsBPos=reshape(Lmax,2,length(Lmax)/2)'; % Insignificant loc max positions
        else
            insCandsB=[];
            insCandsBPos=[];
        end
    end
    clear cands;
    
    % Next
    if c1<tp
        indxStr=sprintf(strg,currentIndex+1);
        eval(['load ',userPath,filesep,'cands',filesep,'cands',indxStr,'.mat;']); % Cands workspace
        candsD=cands; % Structure as is
        Lmax=[cands.Lmax]; 
        allCandsDPos=reshape(Lmax,2,length(Lmax)/2)'; % All local maximum positions
        if ~isempty(cands)
            insCandsD=cands(find([cands.status]==0)); % Cands structure containing only insignificant loc max
            clear Lmax;
            Lmax=[insCandsD.Lmax]; 
            insCandsDPos=reshape(Lmax,2,length(Lmax)/2)'; % Insignificant loc max positions
        else
            insCandsD=[];
            insCandsDPos=[];
        end
    end
    clear cands;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Images
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Load images
    if c1>1
        imgB=imreadnd2(outFileList(c1-1,:),xmin,xmax);
        % Prepare image
        imgB=fsmPrepPrepareImage(imgB,factors(c1-1),W,sigma);
    end
    
    if c1<tp
        imgD=imreadnd2(outFileList(c1+1,:),xmin,xmax);
        % Prepare image
        imgD=fsmPrepPrepareImage(imgD,factors(c1+1),W,sigma);
    end
    
    % For the auxiliary function for the selection of local maxima, the following is needed
    if exist('imgB', 'var')
        % Calculate standard deviation of the non-black pixels as statistical value for comparison
        %nonBlackPixels=imgB(find(imgB>0));
        %stdImg=std(nonBlackPixels(:));
        stdImg = std(imgB(imgB > 0));
        
        % Pad img with intensities equal to the mean img intensity (to avoid border effects)
        maskR=5;
        aImgB=mean(imgB(:))*ones(2*maskR+size(imgB));
        aImgB(maskR+1:end-maskR,maskR+1:end-maskR)=imgB;
    end
    
    if exist('imgD', 'var')
        % Calculate standard deviation of the non-black pixels as statistical value for comparison
        %nonBlackPixels=imgD(find(imgD>0));
        %stdImg=std(nonBlackPixels(:));
        stdImg = std(imgD(imgD > 0));
        
        % Pad img with intensities equal to the mean img intensity (to avoid border effects)
        maskR=5;
        aImgD=mean(imgD(:))*ones(2*maskR+size(imgD));
        aImgD(maskR+1:end-maskR,maskR+1:end-maskR)=imgD;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Fill speckleArray for current frame (c1)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    for c2=1:tm                                 % Cycle through templates
        for c3=1:sp                              % Cycle through speckles (per template and timepoint)
            
            % Reset resp (which accounts for the status 'g' of a speckle)
            resp=0;
            
            % Calculate counter
            counter=c3+(c3-1)*(tp-1)+(c2-1)*(sp*tp)+(c1-1);
            
            % Read indexes
            b=table(counter,1);
            s=table(counter,2);
            d=table(counter,3);
            e=table(counter,4);
           
            if s~=-1
                
                %
                % CASE 1: SPECKLE HAVING NO PRE-BIRTH OR POST-DEATH SPECKLE ASSOCIATED
                %
                if b==0 && d==0
                    speckleArray.timepoint(s)    = uint16(c1);
                    
                    % Check if the speckle is reported in gapList
                    coords(1,1:2)=[cM(c3,(c1-1)*2+1,c2) cM(c3,(c1-1)*2+2,c2)];
                    resp=checkGapList(gapList,coords);
                    speckleArray.spPos(s,1:2)    = uint16(coords'); 
                                        
                    [intensity,background,bk1,bk2,bk3,deltaI,deltaICrit,sigmaMax,sigmaMin,speckleType]=SSpeckleInfo(candsS,allCandsPos,coords);
                    speckleArray.intensity(s)    = intensity;
                    speckleArray.background(s)   = background;
                    speckleArray.bgPos1(s,1:2)   = uint16(bk1');
                    speckleArray.bgPos2(s,1:2)   = uint16(bk2');
                    speckleArray.bgPos3(s,1:2)   = uint16(bk3');
                    speckleArray.deltaI(s)       = deltaI;
                    speckleArray.deltaICrit(s)   = deltaICrit;
                    speckleArray.sigmaSp(s)      = sigmaMax;
                    speckleArray.sigmaBg(s)      = sigmaMin;
                    if resp==1
                        % This speckle is a closed gap
                        speckleArray.status(s)   = 'g';
                    else
                        % This is an actual speckle
                        speckleArray.status(s)   = 's';
                    end
                    if e==1 % This is a speckle in the last frame of the movie
                        speckleArray.status(s)   = 'f';
                    elseif e==2 % This is a speckle in the first frame of the movie
                        speckleArray.status(s)   = 'l';
                    else
                        % Normal speckle
                    end
                    speckleArray.lmEvent(s)      = false;
                    speckleArray.speckleType(s)  = uint8(speckleType);
                    % Add zero score- and activity fields
                    speckleArray.score(s)        = 0;
                    speckleArray.activity(s)     = int8(0);
                end
                
                %
                % CASE 2: SPECKLE HAVING PRE-BIRTH SPECKLE ASSOCIATED
                %
                if b~=0 && d==0
                    
                    % Start with actual speckle
                    speckleArray.timepoint(s)    = uint16(c1);

                    % Check if the speckle is reported in gapList
                    coords(1,1:2)=[cM(c3,(c1-1)*2+1,c2) cM(c3,(c1-1)*2+2,c2)];
                    resp=checkGapList(gapList,coords);
                    speckleArray.spPos(s,1:2)     = uint16(coords'); 
                    
                    [intensity,background,bk1,bk2,bk3,deltaI,deltaICrit,sigmaMax,sigmaMin,speckleType]=SSpeckleInfo(candsS,allCandsPos,coords);
                    speckleArray.intensity(s)    = intensity;
                    speckleArray.background(s)   = background;
                    speckleArray.bgPos1(s,1:2)   = uint16(bk1');
                    speckleArray.bgPos2(s,1:2)   = uint16(bk2');
                    speckleArray.bgPos3(s,1:2)   = uint16(bk3');
                    speckleArray.deltaI(s)       = deltaI;
                    speckleArray.deltaICrit(s)   = deltaICrit;
                    speckleArray.sigmaSp(s)      = sigmaMax;
                    speckleArray.sigmaBg(s)      = sigmaMin;
                    if resp==1
                        % This speckle is a closed gap
                        speckleArray.status(s)   = 'g';
                    else
                        % This is an actual speckle
                        speckleArray.status(s)   = 's';
                    end
                    if e==1 % This is a speckle in the last frame of the movie
                        speckleArray.status(s)   = 'f';
                    elseif e==2 % This is a speckle in the first frame of the movie
                        speckleArray.status(s)   = 'l';
                    else
                        % Normal speckle
                    end
                    speckleArray.lmEvent(s)      = false;
                    speckleArray.speckleType(s)  = uint8(speckleType);
                    % Add zero score- and activity fields
                    speckleArray.score(s)        = 0;
                    speckleArray.activity(s)     = int8(0);
                    
                    % Now with pre-birth non-speckle
                    speckleArray.timepoint(b)    = speckleArray.timepoint(s)-1; % This requires an overloaded operator for uint16
                    [intensity,background,bk1,bk2,bk3,deltaI,deltaICrit,sigmaMax,sigmaMin,lmEvent,spPos,speckleType]=BDSpeckleInfo(coords,candsS,allCandsPos,insCandsB,insCandsBPos,imgB,c1,c2,c3,cM,noiseParams,threshold,aImgB,maskR,stdImg);
                    speckleArray.spPos(b,1:2)    = uint16(spPos');            
                    speckleArray.intensity(b)    = intensity;
                    speckleArray.background(b)   = background;
                    speckleArray.bgPos1(b,1:2)   = uint16(bk1');
                    speckleArray.bgPos2(b,1:2)   = uint16(bk2');
                    speckleArray.bgPos3(b,1:2)   = uint16(bk3');
                    speckleArray.deltaI(b)       = deltaI;
                    speckleArray.deltaICrit(b)   = deltaICrit;
                    speckleArray.sigmaSp(b)      = sigmaMax;
                    speckleArray.sigmaBg(b)      = sigmaMin;
                    speckleArray.status(b)       = 'b';
                    speckleArray.lmEvent(b)      = logical(lmEvent);
                    speckleArray.speckleType(b)  = uint8(speckleType);
                    % Add zero score- and activity fields
                    speckleArray.score(b)        = 0;
                    speckleArray.activity(b)     = int8(0);
                end
                
                %
                % CASE 3: SPECKLE HAVING POST-DEATH SPECKLE ASSOCIATED
                %
                if b==0 && d~=0
                    
                    % Start with actual speckle
                    speckleArray.timepoint(s)    = uint16(c1);

                    % Check if the speckle is reported in gapList
                    coords(1,1:2)=[cM(c3,(c1-1)*2+1,c2) cM(c3,(c1-1)*2+2,c2)];
                    resp=checkGapList(gapList,coords);
                    speckleArray.spPos(s,1:2)    = uint16(coords'); 
                    
                    [intensity,background,bk1,bk2,bk3,deltaI,deltaICrit,sigmaMax,sigmaMin,speckleType]=SSpeckleInfo(candsS,allCandsPos,coords);
                    speckleArray.intensity(s)    = intensity;
                    speckleArray.background(s)   = background;
                    speckleArray.bgPos1(s,1:2)   = uint16(bk1');
                    speckleArray.bgPos2(s,1:2)   = uint16(bk2');
                    speckleArray.bgPos3(s,1:2)   = uint16(bk3');
                    speckleArray.deltaI(s)       = deltaI;
                    speckleArray.deltaICrit(s)   = deltaICrit;
                    speckleArray.sigmaSp(s)      = sigmaMax;
                    speckleArray.sigmaBg(s)      = sigmaMin;
                    if resp==1
                        % This speckle is a closed gap
                        speckleArray.status(s)   = 'g';
                    else
                        % This is an actual speckle
                        speckleArray.status(s)   = 's';
                    end
                    if e==1 % This is a speckle in the last frame of the movie
                        speckleArray.status(s)   = 'f';
                    elseif e==2 % This is a speckle in the first frame of the movie
                        speckleArray.status(s)   = 'l';
                    else
                        % Normal speckle
                    end
                    speckleArray.lmEvent(s)      = false;
                    speckleArray.speckleType(s)  = uint8(speckleType);
                    % Add zero score- and activity fields
                    speckleArray.score(s)        = 0;
                    speckleArray.activity(s)     = int8(0);
                    
                    % Now with post-death non-speckle
                    speckleArray.timepoint(d)    = speckleArray.timepoint(s)+1; % This requires an overloaded operator for uint16
                    [intensity,background,bk1,bk2,bk3,deltaI,deltaICrit,sigmaMax,sigmaMin,lmEvent,spPos,speckleType]=BDSpeckleInfo(coords,candsS,allCandsPos,insCandsD,insCandsDPos,imgD,c1,c2,c3,cM,noiseParams,threshold,aImgD,maskR,stdImg);
                    speckleArray.spPos(d,1:2)    = uint16(spPos');         
                    speckleArray.intensity(d)    = intensity;
                    speckleArray.background(d)   = background;
                    speckleArray.bgPos1(d,1:2)   = uint16(bk1');
                    speckleArray.bgPos2(d,1:2)   = uint16(bk2');
                    speckleArray.bgPos3(d,1:2)   = uint16(bk3');
                    speckleArray.deltaI(d)       = deltaI;
                    speckleArray.deltaICrit(d)   = deltaICrit;
                    speckleArray.sigmaSp(d)      = sigmaMax;
                    speckleArray.sigmaBg(d)      = sigmaMin;
                    speckleArray.status(d)       = 'd';
                    speckleArray.lmEvent(d)      = logical(lmEvent);			  
                    speckleArray.speckleType(d)  = uint8(speckleType);
                    % Add zero score- and activity fields
                    speckleArray.score(d)        = 0;
                    speckleArray.activity(d)     = int8(0);
                    
                end
                
                %
                % CASE 4: SPECKLE HAVING BOTH PRE-BIRTH OR POST-DEATH SPECKLE ASSOCIATED
                %
                if b~=0 && d~=0                             % THIS IS A GHOST SPECKLE!
                    
                    % Start with actual speckle
                    speckleArray.timepoint(s)    = uint16(c1);

                    % Check if the speckle is reported in gapList
                    coords(1,1:2)=[cM(c3,(c1-1)*2+1,c2) cM(c3,(c1-1)*2+2,c2)];
                    resp=checkGapList(gapList,coords);
                    speckleArray.spPos(s,1:2)    = uint16(coords'); 
                    
                    [intensity,background,bk1,bk2,bk3,deltaI,deltaICrit,sigmaMax,sigmaMin,speckleType]=SSpeckleInfo(candsS,allCandsPos,coords);
                    speckleArray.intensity(s)    = intensity;
                    speckleArray.background(s)   = background;
                    speckleArray.bgPos1(s,1:2)   = uint16(bk1');
                    speckleArray.bgPos2(s,1:2)   = uint16(bk2');
                    speckleArray.bgPos3(s,1:2)   = uint16(bk3');
                    speckleArray.deltaI(s)       = deltaI;
                    speckleArray.deltaICrit(s)   = deltaICrit;
                    speckleArray.sigmaSp(s)      = sigmaMax;
                    speckleArray.sigmaBg(s)      = sigmaMin;
                    if resp==1
                        % This speckle is a closed gap
                        speckleArray.status(s)   = 'g';
                    else
                        % This is an actual speckle
                        speckleArray.status(s)   = 's';
                    end
                    if e==1 % This is a speckle in the last frame of the movie
                        speckleArray.status(s)   = 'f';
                    elseif e==2 % This is a speckle in the first frame of the movie
                        speckleArray.status(s)   = 'l';
                    else
                        % Normal speckle
                    end
                    speckleArray.lmEvent(s)      = false;
                    speckleArray.speckleType(s)  = uint8(speckleType);
                    % Add zero score- and activity fields
                    speckleArray.score(s)        = 0;
                    speckleArray.activity(s)     = int8(0);
                    
                    % Now with pre-birth non-speckle
                    speckleArray.timepoint(b)    = speckleArray.timepoint(s)-1; % This requires an overloaded operator for uint16
                    [intensity,background,bk1,bk2,bk3,deltaI,deltaICrit,sigmaMax,sigmaMin,lmEvent,spPos,speckleType]=BDSpeckleInfo(coords,candsS,allCandsPos,insCandsB,insCandsBPos,imgB,c1,c2,c3,cM,noiseParams,threshold,aImgB,maskR,stdImg);
                    speckleArray.spPos(b,1:2)    = uint16(spPos');            
                    speckleArray.intensity(b)    = intensity;
                    speckleArray.background(b)   = background;
                    speckleArray.bgPos1(b,1:2)   = uint16(bk1');
                    speckleArray.bgPos2(b,1:2)   = uint16(bk2');
                    speckleArray.bgPos3(b,1:2)   = uint16(bk3');
                    speckleArray.deltaI(b)       = deltaI;
                    speckleArray.deltaICrit(b)   = deltaICrit;
                    speckleArray.sigmaSp(b)      = sigmaMax;
                    speckleArray.sigmaBg(b)      = sigmaMin;
                    speckleArray.status(b)       = 'b';
                    speckleArray.lmEvent(b)      = logical(lmEvent);
                    speckleArray.speckleType(b)  = uint8(speckleType);
                    % Add zero score- and activity fields
                    speckleArray.score(b)        = 0;
                    speckleArray.activity(b)     = int8(0);
                    
                    % Now with post-death non-speckle
                    speckleArray.timepoint(d)    = speckleArray.timepoint(s)+1; % This requires an overloaded operator for uint16
                    [intensity,background,bk1,bk2,bk3,deltaI,deltaICrit,sigmaMax,sigmaMin,lmEvent,spPos,speckleType]=BDSpeckleInfo(coords,candsS,allCandsPos,insCandsD,insCandsDPos,imgD,c1,c2,c3,cM,noiseParams,threshold,aImgD,maskR,stdImg);
                    speckleArray.spPos(d,1:2)    = uint16(spPos');
                    speckleArray.intensity(d)    = intensity;
                    speckleArray.background(d)   = background;
                    speckleArray.bgPos1(d,1:2)   = uint16(bk1');
                    speckleArray.bgPos2(d,1:2)   = uint16(bk2');
                    speckleArray.bgPos3(d,1:2)   = uint16(bk3');             
                    speckleArray.deltaI(d)       = deltaI;
                    speckleArray.deltaICrit(d)   = deltaICrit;
                    speckleArray.sigmaSp(d)      = sigmaMax;
                    speckleArray.sigmaBg(d)      = sigmaMin;
                    speckleArray.status(d)       = 'd';
                    speckleArray.lmEvent(d)      = logical(lmEvent);
                    speckleArray.speckleType(d)  = uint8(speckleType);
                    % Add zero score- and activity fields
                    speckleArray.score(d)        = 0;
                    speckleArray.activity(d)     = int8(0);
                end  
            end
            
            % Update waitbar if needed
            if step>1
                if mod((c1*c2*c3),step)==1
                    waitbar((c1*c2*c3)/tot,h);
                end
            else
                waitbar((c1*c2*c3)/tot,h);
            end        
            
        end
    end
end

% Close waitbar
close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [intensity,background,bk1,bk2,bk3,deltaI,deltaICrit,sigmaMax,sigmaMin,speckleType]=SSpeckleInfo(candsS,allCandsPos,lmPos)

% Look for 'coords' in the cands structure
t=lmPos(1,1)==allCandsPos(:,1);
u=lmPos(1,2)==allCandsPos(:,2);
v=find(t & u);

% Check for correct selection of the loc max
if isempty(v) || length(v)>1
    error('Loc max not univocally found in Delaunay Triangulation results');
end

% This is a speckle of class 's' -> all information has been already stored in delaunay###.mat
intensity=candsS(v).ILmax;
background=candsS(v).IBkg;
bk1=candsS(v).Bkg1;
bk2=candsS(v).Bkg2;
bk3=candsS(v).Bkg3;
deltaI=candsS(v).deltaI;
deltaICrit=candsS(v).deltaICrit;
sigmaMax=candsS(v).sigmaLmax;
sigmaMin=candsS(v).sigmaBkg;
speckleType=candsS(v).speckleType;


function [intensity,background,bk1,bk2,bk3,deltaI,deltaICrit,sigmaMax,sigmaMin,lmEvent,spPos,speckleType]=BDSpeckleInfo(lmPos,candsS,allCandsPos,candsEv,candsEvPos,img,timepoint,c2,c3,cM,noiseParams,threshold,aImg,maskR,stdImg)

% Look for 'coords' in the cands structure (still for the 's' speckle)
t=lmPos(1,1)==allCandsPos(:,1);
u=lmPos(1,2)==allCandsPos(:,2);
v=find(t & u);

% Check for correct selection of the loc max
if isempty(v) || length(v)>1
    error('Loc max not univocally found in Delaunay Triangulation results');
end

%
% NOW LOOK FOR A WEAK LOCAL MAXIMUM FOR THE 'b' AND 'd' SPECKLE; OTHERWISE CALCULATE
%     NOISE VALUES AS USUAL
%

if isempty(lmPos) || isempty(candsEvPos)
    D=[];
else
    D=createDistanceMatrix(lmPos,candsEvPos);
end

if ~isempty(D)
    
    % Row
    E=zeros(size(D));
    for i=1:size(D,1)
        t=D(i,:); 
        t=t==min(t); 
        E(i,:)=t;
    end
    
    % Column
    F=zeros(size(D));
    for i=1:size(D,2)
        t=D(:,i); 
        t=t==min(t); 
        F(:,i)=t;
    end
    
    % Thresholding
    H=D<=threshold;
    
    % Resulting selected distance matrix
    G=(E & F & H);
    
    % Now find matching pairs: y = position in lmPos (must be 1), x = position in evPos
    [y x]=find(G);
    
else
    
    y=[]; x=[];
    
end

if ~isempty(y) && ~isempty(x)
    if y~=1
        error('fsmBuildSaveSpeckleArray::BDSpeckleInfo - Too many source speckles!');
    else
        if length(x)>1
            % Arbitrarily take the one with smaller deltaI
            intensity = img(sub2ind(size(img), candsEvPos(x(:), 1), ...
                candsEvPos(x(:), 2)))'; 
            background = cell2mat({candsEv(x(:)).IBkg});
            deltaI = intensity - background;
            index=find(deltaI==min(deltaI)); 
            if length(index)>1
                disp('Arbitrarily chosen one local maximum as a ''b'' or ''d'' event.');
                index=index(1);
            end
            x=x(index);
        end
        % Read intensities at the coordinates of the local maximum found
        intensity=candsEv(x).ILmax;
        background=candsEv(x).IBkg;
        bk1=candsEv(x).Bkg1;
        bk2=candsEv(x).Bkg2;
        bk3=candsEv(x).Bkg3;
        deltaI=candsEv(x).deltaI;
        deltaICrit=candsEv(x).deltaICrit;
        sigmaMax=candsEv(x).sigmaLmax;
        sigmaMin=candsEv(x).sigmaBkg;
        lmEvent=1; % Insignificant local maximum found
        speckleType=candsEv(x).speckleType;
    end
else % No matching discarded local maximum found
    speckleType=0; % This is not even a weak local maximum
    % Read intensities at the coordinates of speckle 's'
    intensity=img(lmPos(1),lmPos(2));
    % Before using the data of the speckle 's' to read the corresponding data for 'b' or 'd'
    %    one has to check whether the loc max has been validated using Delaunay or the auxiliary
    %    function
    if candsS(v).Bkg1(1)==-1 % Auxiliary function
        disp('Using auxiliary function to recover speckle info for non-existant loc max (failed Delaunay triangulation)');
        [background,deltaI,k,sigmaDiff,sigmaMax,sigmaMin,status,A]=fsmPrepTestSpeckleSignifAux(aImg,img,[lmPos(1) lmPos(2)],maskR,stdImg,0,noiseParams);
        bk1=candsS(v).Bkg1;
        bk2=candsS(v).Bkg2;
        bk3=candsS(v).Bkg3;
        lmEvent=0; % No insignificant local maximum found
        deltaICrit=k*sigmaDiff;	
    else % Delaunay
        background=mean([img(candsS(v).Bkg1(1),candsS(v).Bkg1(2)) img(candsS(v).Bkg2(1),candsS(v).Bkg2(2)) img(candsS(v).Bkg3(1),candsS(v).Bkg3(2))]);
        bk1=candsS(v).Bkg1;
        bk2=candsS(v).Bkg2;
        bk3=candsS(v).Bkg3;
        lmEvent=0; % No insignificant local maximum found
        [Imin,deltaI,k,sigmaDiff,sigmaMax,sigmaMin,status,A]=fsmPrepTestSpeckleSignif([lmPos(1) lmPos(2)],intensity,background,noiseParams(1),noiseParams(2),noiseParams(3),noiseParams(4),0);
        deltaICrit=k*sigmaDiff;
    end
    
end

% Return also the position of the recovered local maximum
if ~isempty(x)
    spPos=candsEvPos(x,1:2);
else
    spPos=lmPos; % If nothing has been found, use the coordinates of the 's' speckle
end

%
%
%
%
%

function resp=checkGapList(gapList,coords)

% Check
if size(coords,1)>1 || size(coords,1)==0
    error('Either no speckle or more than one speckle passed to fsmBuildSaveSpeckleArray::checkGapList');
end
t=coords(1,1)==gapList(:,1);
u=coords(1,2)==gapList(:,2);
v=find(t & u);
if isempty(v)
    resp=0;
end
if length(v)==1
    resp=1;
end
if length(v)>1
    error('fsmBuildSaveSpeckleArray: speckle found more than once in gapList!');
end

