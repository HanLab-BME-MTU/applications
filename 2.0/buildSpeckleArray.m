function speckleArray=buildSpeckleArray(cM,noiseParams,threshold,gapList,cands,img,varargin)
% buildSpeckleArray fills the speckleArray structure containing all speckle information extracted from the analysis of an FSM stack
%
% SYNOPSIS   speckleArray=buildSpeckleArray(cM,noiseParams,threshold)
%
% INPUT      cM            : magic position matrix (mpm) as returned by the tracker module
%            noiseParams   : noise parameters for speckle significance test
%            threshold     : radius of the region considered when looking for matching speckles
%
%
% OUTPUT     speckleArray  : structure containing all speckle information extracted from an image stack

% Sebastien Besson, June 2011
% Adapted from fsmBuildSaveSpeckleArray

ip =inputParser;
ip.addRequired('cM',@isnumeric);
ip.addRequired('noiseParams',@isnumeric);
ip.addRequired('threshold',@isscalar);
ip.addRequired('gapList',@iscell);
ip.addRequired('cands',@iscell);
ip.addRequired('img',@iscell);
ip.addParamValue('waitbar',[],@ishandle);
ip.parse(cM,noiseParams,threshold,gapList,cands,img,varargin{:});


nSpeckles=size(cM,1);   % Number of speckle lines
nFrames=size(cM,2)/2; % Number of time points
nTemplates=size(cM,3);   % Number of subwindows

% Initialize waitbar 
if ~isempty(ip.Results.waitbar)
    wtBar=ip.Results.waitbar;
    waitbar(0,wtBar,'Indexing speckles...');
elseif feature('ShowFigureWindows')
    wtBar=waitbar(0,'Indexing speckles...');
end

% Index speckles and non-speckles
x=1:2:2*nFrames-1;

% Allocate memory for 'table'
table=zeros(nSpeckles*nFrames,4);

pos4=0;                               % This marks a speckle still existing in the last time point
cMPos=0;                              % Position in cM
sAPos=0;                              % Position in speckleArray
for iSpeckle=1:nSpeckles
    for iFrame=1:nFrames
        cMPos=cMPos+1;                  % Keep info of the position in cM  
        currentSpeckle=cM(iSpeckle,x(iFrame));    % Read x coordinate of current speckle
        
        if currentSpeckle==0
            pos1=-1;
            pos2=-1;
            pos3=-1;
        else
            % We are at the beginning of the movie
            if iFrame==1
                pos4=1;                  % We have to mark this because we have a speckle already existing at the first time point
                if cM(iSpeckle,x(iFrame+1))~=0      % Also next one is speckle
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
                
                %
                % We are at the end of the movie
                %
            elseif iFrame==nFrames
                pos4=2;                  % We have to cut here because we are at the last time point with a still existing speckle
                if cM(iSpeckle,x(iFrame-1))==0      % Last one was a non-speckle
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
                
                %
                % We are in the central time-points
                %
                
            else
                
                if cM(iSpeckle,x(iFrame-1))~=0 && cM(iSpeckle,x(iFrame+1))~=0   % Last and next ones are speckle
                    pos1=0;
                    sAPos=sAPos+1;          % This is the speckle
                    pos2=sAPos;
                    pos3=0;
                end
                
                if cM(iSpeckle,x(iFrame-1))==0 && cM(iSpeckle,x(iFrame+1))~=0   % Last one was not a speckle; the next is
                    sAPos=sAPos+1;
                    pos1=sAPos;             % This is the "pre-birth" non-speckle
                    sAPos=sAPos+1;
                    pos2=sAPos;             % This is the speckle
                    pos3=0;
                end
                
                if cM(iSpeckle,x(iFrame-1))~=0 && cM(iSpeckle,x(iFrame+1))==0   % Last one was a speckle; the next is not
                    pos1=0;
                    sAPos=sAPos+1;          % This is the speckle
                    pos2=sAPos;
                    sAPos=sAPos+1;          % This is the "post-death" non-speckle
                    pos3=sAPos;
                end
                
                if cM(iSpeckle,x(iFrame-1))==0 && cM(iSpeckle,x(iFrame+1))==0   % Both are not speckles
                    sAPos=sAPos+1;
                    pos1=sAPos;             % This is the "pre-birth" non-speckle
                    sAPos=sAPos+1;
                    pos2=sAPos;             % This is the speckle
                    sAPos=sAPos+1;
                    pos3=sAPos;             % This is the "post-death" non-speckle
                end
                
            end
        end
        % Write correspondence
        table(cMPos,1:4)=[pos1 pos2 pos3 pos4];
        % Reset pos4;
        pos4=0;
                
    end
end


% Initialize empty speckle array
lengthSpeckleArray=max(table(:));
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
% For vectorization
% speckleArray2=speckleArray;

% Initialize waitbar 
if ishandle(wtBar), 
    waitbar(0,wtBar,'Populating speckleArray structure...'); 
end
nTot = nSpeckles*nFrames*nTemplates;
for iFrame=1:nFrames                                    % Cycle through timepoints
    
    % Current
    candsS=cands{iFrame}; % Structure as is

    % Previous
    if iFrame>1, 
        insCandsB=cands{iFrame-1}([cands{iFrame-1}.status]==0); 
        imgB =img{iFrame-1};
    end
    
    % Next
    
    if iFrame<nFrames
        insCandsD=cands{iFrame+1}([cands{iFrame+1}.status]==0);
        imgD =img{iFrame+1};
    end
    
%     %%%% Beginning of vectorization
%     tic
%     [iTemplate,iSpeckle] =meshgrid(1:nTemplates,1:nSpeckles);
%     counter=iSpeckle+(iSpeckle-1)*(nFrames-1)+(iTemplate-1)*(nSpeckles*nFrames)+(iFrame-1);
%     % Read indexes
%     b=table(counter,1);
%     s=table(counter,2);
%     d=table(counter,3);
%     e=table(counter,4);
%     
%     % Valid speckles
%     validIndx=(s~=-1);
%     validS=s(validIndx);
%     speckleArray2.timepoint(validS)=uint16(iFrame);
%     coords2(:,1) =arrayfun(@(x,y) cM(x,(iFrame-1)*2+1,y),...
%         iSpeckle(validIndx),iTemplate(validIndx));
%     coords2(:,2) =arrayfun(@(x,y) cM(x,(iFrame-1)*2+2,y),...
%         iSpeckle(validIndx),iTemplate(validIndx));
% %     coords2 = speckleArray2.spPos(validS,:);
%     speckleArray2.spPos(validS,:) =uint16(coords2);
%     
%     % FInd associated speckles candidates
%     allCandsSPos = vertcat(candsS.Lmax);
%     [~,icoords,icandsPos]=intersect(coords2,allCandsSPos,'rows');
%     [sortedIcoords,indx] = sort(icoords);
%     if ~isequal(sortedIcoords',1:size(coords2,1)),
%         error('Loc max not univocally found in Delaunay Triangulation results');
%     end
%     
%     % Read various speckle data
%     speckleArray2.intensity(validS) = [candsS(icandsPos(indx)).ILmax];
%     speckleArray2.background(validS) = [candsS(icandsPos(indx)).IBkg];
%     speckleArray2.bgPos1(validS,:) = uint16(vertcat(candsS(icandsPos(indx)).Bkg1));
%     speckleArray2.bgPos2(validS,:) = uint16(vertcat(candsS(icandsPos(indx)).Bkg2));
%     speckleArray2.bgPos3(validS,:) = uint16(vertcat(candsS(icandsPos(indx)).Bkg3));
%     speckleArray2.deltaI(validS) = [candsS(icandsPos(indx)).deltaI];
%     speckleArray2.deltaICrit(validS) = [candsS(icandsPos(indx)).deltaICrit];
%     speckleArray2.sigmaSp(validS) = [candsS(icandsPos(indx)).sigmaLmax];
%     speckleArray2.sigmaBg(validS) = [candsS(icandsPos(indx)).sigmaBkg];
%     speckleArray2.speckleType(validS) = uint8([candsS(icandsPos(indx)).speckleType]);
%     speckleArray2.lmEvent(validS) = false;
%     % Set the status to 's' by default
%     speckleArray2.status(validS)='s';
%     speckleArray2.status(s(e==1))='f';
%     speckleArray2.status(s(e==2))='l';
%     [~,icoords]=intersect(coords2,gapList{iFrame},'rows');
%     speckleArray2.status(validS(icoords))='g';
%     
%     
%     % Birth events
%     validB=b(validIndx);
%     birthEvents = (validB~=0);
%     if any(birthEvents)
%         speckleArray2.timepoint(validB(birthEvents))=...
%             speckleArray2.timepoint(validS(birthEvents))-1;
%         speckleArray2.status(validB(birthEvents))='b';
%         
%         [intensity,background,bk1,bk2,bk3,deltaI,deltaICrit,sigmaMax,sigmaMin,lmEvent,spPos,speckleType]=...
%             associateSpeckles(coords2(birthEvents,:),insCandsB,imgB,noiseParams,threshold);
% 
% %         insCandsBPos=vertcat(insCandsB.Lmax);
% %         D=KDTreeBallQuery(insCandsBPos,coords2(birthEvents,:),threshold);
%     end
%     
%     % Death events
%     potentialD=d(validIndx);
%     deathEvents = (potentialD~=0);
%     validD=potentialD(deathEvents);
%     if any(deathEvents)
%         speckleArray2.timepoint(validD)=...
%             speckleArray2.timepoint(validS(deathEvents))+1;
%         speckleArray2.status(validD)='d';
% 
%         %
%         % NOW LOOK FOR A WEAK LOCAL MAXIMUM FOR THE 'b' AND 'd' SPECKLE; OTHERWISE CALCULATE
%         %     NOISE VALUES AS USUAL
%         %
%         candsEvPos=vertcat(insCandsD.Lmax);
%         [D dist]=KDTreeBallQuery(candsEvPos,coords2(deathEvents,:),threshold);
%         unassociatedD = cellfun(@isempty,D);
%         speckleArray2.lmEvent(validD(unassociatedD))=0; % This is not even a weak local maximum
%         speckleArray2.lmEvent(validD(~unassociatedD))=1;
%         % Read intensities at the coordinates of speckle 's'
% %         speckleArray2.speckleType(deathEvents&unassociatedD) =imgD(lmPos(1),lmPos(2));
% %         % Before using the data of the speckle 's' to read the corresponding data for 'b' or 'd'
% %         %    one has to check whether the loc max has been validated using Delaunay or the auxiliary
% %         %    function
% %         
% %         
% %         bk1=cands(v).Bkg1;
% %         bk2=cands(v).Bkg2;
% %         bk3=cands(v).Bkg3;
% %         % SB: some background points can be outside the cell mask where
% %         % image is null. Call nonzeros to compute mean over non zeros
% %         % pixels.
% %         background=mean(nonzeros([img(bk1(1),bk1(2)) img(bk2(1),bk2(2)) ...
% %             img(bk3(1),bk3(2))]));
% %         lmEvent=0; % No insignificant local maximum found
% %         [Imin,deltaI,k,sigmaDiff,sigmaMax,sigmaMin,status]=testSpeckleSignificance(intensity,background,noiseParams(1),noiseParams(2),noiseParams(3),noiseParams(4));
% %         deltaICrit=k*sigmaDiff;
% % %         [intensity,background,bk1,bk2,bk3,deltaI,deltaICrit,sigmaMax,sigmaMin,lmEvent,spPos,speckleType]=...
% %             associateSpeckles(coords2(deathEvents,:),insCandsD,imgD,noiseParams,threshold);
%     end
% %         
% %         coordsD = coords2(deathEvents,:);
% %         insCandsDPos=vertcat(insCandsD.Lmax);
% %         D=KDTreeBallQuery(insCandsDPos,coordsD,...
% %             threshold*ones(sum(deathEvents),1));
% %         
% %         % Death events not-matching any insignificant local maxima
% %         umDEvents =cellfun(@isempty,D);
% %         validUmD=validD(umDEvents);
% %         speckleArray2.intensity(validUmD) = ...
% %             arrayfun(@(x)imgD(coordsD(x,1),coordsD(x,2)),find(umDEvents));
% %         speckleArray2.speckleType(validUmD) = 0;
% %         
% %         [~,icoords,icandsPos]=intersect(coordsD(umDEvents,:),allCandsSPos,'rows');
% %         speckleArray2.bgPos1(validUmD,:) = vertcat(candsS(icandsPos).Bkg1);
% %         speckleArray2.bgPos2(validUmD,:) = vertcat(candsS(icandsPos).Bkg2);
% %         speckleArray2.bgPos3(validUmD,:) = vertcat(candsS(icandsPos).Bkg3);
% %         speckleArray2.lmEvent(validUmD) = 0; % No insignificant local maximum found
% %         % SB: some background points can be outside the cell mask where
% %         % image is null. Call nonzeros to compute mean over non zeros
% %         % pixels.
% %         %         speckleArray2.background(validD(umDEvents)) = arrayfun(@(x)...
% %         %             mean(nonzeros([img(bk1(1),bk1(2)) img(bk2(1),bk2(2)) ...
% %         % %             img(bk3(1),bk3(2))]));
% %         speckleArray2.background(validUmD) = arrayfun(@(x)...
% %             mean(nonzeros([imgD(candsS(x).Bkg1(1),candsS(x).Bkg1(2))...
% %             imgD(candsS(x).Bkg2(1),candsS(x).Bkg2(2))...
% %             imgD(candsS(x).Bkg3(1),candsS(x).Bkg3(2))])),icandsPos);
% %         %         [Imin,deltaI,k,sigmaDiff,sigmaMax,sigmaMin,status]=testSpeckleSignificance(intensity,background,noiseParams(1),noiseParams(2),noiseParams(3),noiseParams(4));
% %         %         deltaICrit=k*sigmaDiff;
% %         
% %         % Death events matching insignificant local maxima
% %         validMD=validD(~umDEvents);
% %         speckleArray2.lmEvent(validMD) = 1;
% %     end
% %     toc
%     %%%% End of vectorization

 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Fill speckleArray for current frame (iFrame)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for iTemplate=1:nTemplates                   % Cycle through templates
        for iSpeckle=1:nSpeckles  % Cycle through speckles (per template and timepoint)
           
            % Calculate counter
            counter=iSpeckle+(iSpeckle-1)*(nFrames-1)+(iTemplate-1)*(nSpeckles*nFrames)+(iFrame-1);
            
            % Read indexes
            b=table(counter,1);
            s=table(counter,2);
            d=table(counter,3);
            e=table(counter,4);
           
            if s~=-1
                 
                % Start with actual speckle
                speckleArray.timepoint(s)    = uint16(iFrame);
                
                coords(1,1:2)=[cM(iSpeckle,(iFrame-1)*2+1,iTemplate) cM(iSpeckle,(iFrame-1)*2+2,iTemplate)];
                speckleArray.spPos(s,1:2)     = uint16(coords');
                
                [intensity,background,bk1,bk2,bk3,deltaI,deltaICrit,sigmaMax,sigmaMin,speckleType]=SSpeckleInfo(candsS,coords);
                speckleArray.intensity(s)    = intensity;
                speckleArray.background(s)   = background;
                speckleArray.bgPos1(s,1:2)   = uint16(bk1');
                speckleArray.bgPos2(s,1:2)   = uint16(bk2');
                speckleArray.bgPos3(s,1:2)   = uint16(bk3');
                speckleArray.deltaI(s)       = deltaI;
                speckleArray.deltaICrit(s)   = deltaICrit;
                speckleArray.sigmaSp(s)      = sigmaMax;
                speckleArray.sigmaBg(s)      = sigmaMin;
                if ismember(coords,gapList{iFrame},'rows')
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
                end
                speckleArray.lmEvent(s)      = false;
                speckleArray.speckleType(s)  = uint8(speckleType);
                
                % SPECKLE HAVING PRE-BIRTH SPECKLE ASSOCIATED
                if b~=0   
                    
                    % Now with pre-birth non-speckle
                    speckleArray.timepoint(b)    = speckleArray.timepoint(s)-1; % This requires an overloaded operator for uint16
                    [intensity,background,bk1,bk2,bk3,deltaI,deltaICrit,sigmaMax,sigmaMin,lmEvent,spPos,speckleType]=BDSpeckleInfo(coords,candsS,insCandsB,imgB,noiseParams,threshold);
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
                end
                
                % CASE 3: SPECKLE HAVING POST-DEATH SPECKLE ASSOCIATED
                if d
                    % Now with post-death non-speckle
                    speckleArray.timepoint(d)    = speckleArray.timepoint(s)+1; % This requires an overloaded operator for uint16
                    [intensity,background,bk1,bk2,bk3,deltaI,deltaICrit,sigmaMax,sigmaMin,lmEvent,spPos,speckleType]=BDSpeckleInfo(coords,candsS,insCandsD,imgD,noiseParams,threshold);
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

                    
                end
%                 
            end
            
            
            % Update waitbar if needed
            linInd = sub2ind([nSpeckles nFrames nTemplates],iSpeckle,iFrame,iTemplate);
            if mod(linInd,round(nTot/20))==1 && ishandle(wtBar)
                ti=linInd/nTot;
                waitbar(ti,wtBar);
            end
            
        end
    end
    
%     toc
%     isequal(speckleArray,speckleArray2)
end

% Close waitbar if not-delegated
if isempty(ip.Results.waitbar) && ishandle(wtBar), 
    close(wtBar); 
end

function [intensity,background,bk1,bk2,bk3,deltaI,deltaICrit,sigmaMax,sigmaMin,speckleType]=SSpeckleInfo(cands,lmPos)

allCandsPos=vertcat(cands.Lmax);
% Look for 'coords' in the cands structure
t=lmPos(1,1)==allCandsPos(:,1);
u=lmPos(1,2)==allCandsPos(:,2);
v=find(t & u);

% Check for correct selection of the loc max
if isempty(v) || length(v)>1
    error('Loc max not univocally found in Delaunay Triangulation results');
end

% This is a speckle of class 's' -> all information has been already stored in delaunay###.mat
intensity=cands(v).ILmax;
background=cands(v).IBkg;
bk1=cands(v).Bkg1;
bk2=cands(v).Bkg2;
bk3=cands(v).Bkg3;
deltaI=cands(v).deltaI;
deltaICrit=cands(v).deltaICrit;
sigmaMax=cands(v).sigmaLmax;
sigmaMin=cands(v).sigmaBkg;
speckleType=cands(v).speckleType;


function [intensity,background,bk1,bk2,bk3,deltaI,deltaICrit,sigmaMax,sigmaMin,lmEvent,spPos,speckleType]=BDSpeckleInfo(lmPos,cands,candsEv,img,noiseParams,threshold)

allCandsPos=vertcat(cands.Lmax);
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
candsEvPos=vertcat(candsEv.Lmax);
%SB: quick
if isempty(candsEvPos), x=[];
else
    D=KDTreeBallQuery(candsEvPos,lmPos,threshold);
    
    % To be compatibile with the fsmBuildSaveSpeckleArray
    % index=(dist{1}==min(dist{1}));
    % x=D{1}(index);
    x=D{1};
end

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

if ~isempty(x)
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
else % No matching discarded local maximum found
    speckleType=0; % This is not even a weak local maximum
    % Read intensities at the coordinates of speckle 's'
    intensity=img(lmPos(1),lmPos(2));
    % Before using the data of the speckle 's' to read the corresponding data for 'b' or 'd'
    %    one has to check whether the loc max has been validated using Delaunay or the auxiliary
    %    function
    
    
    bk1=cands(v).Bkg1;
    bk2=cands(v).Bkg2;
    bk3=cands(v).Bkg3;
    % SB: some background points can be outside the cell mask where
    % image is null. Call nonzeros to compute mean over non zeros
    % pixels.
    background=mean(nonzeros([img(bk1(1),bk1(2)) img(bk2(1),bk2(2)) ...
        img(bk3(1),bk3(2))]));
    lmEvent=0; % No insignificant local maximum found
    [Imin,deltaI,k,sigmaDiff,sigmaMax,sigmaMin,status]=testSpeckleSignificance(intensity,background,noiseParams(1),noiseParams(2),noiseParams(3),noiseParams(4));
    deltaICrit=k*sigmaDiff;
end

% Return also the position of the recovered local maximum
if ~isempty(x)
    spPos=candsEvPos(x,1:2);
else
    spPos=lmPos; % If nothing has been found, use the coordinates of the 's' speckle
end

