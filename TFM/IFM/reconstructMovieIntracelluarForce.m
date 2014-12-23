function reconstructMovieIntracelluarForce(movieData,varargin)
% a code for intracellular force reconstruction
% 
%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous flow tracking processes                                                                     
iFlowProc = movieData.getProcessIndex('FlowTrackingProcess',1,0);
p = movieData.getProcess(iFlowProc).funParams_;

% iProc = movieData.getProcessIndex('FlowAnalysisProcess',1,0);
% %If the process doesn't exist, create it
% if isempty(iProc)
%     iProc = numel(movieData.processes_)+1;
%     movieData.addProcess(FlowAnalysisProcess(movieData,...
%         movieData.outputDirectory_));                                                                                                 
% end
% 
% flowAnProc = movieData.processes_{iProc};
% %Parse input, store in parameter structure
% p = parseProcessParams(flowAnProc,paramsIn);
% if strcmp(p.FlowProcess,'SpeckleTrackingProcess')
%     p.FlowProcess = 'FlowTrackingProcess';
% end

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name','Intracellular force reconstruction ...');
end

% Reading various constants
nFrames = movieData.nFrames_;

% Test the presence and output validity of the speckle tracking process
iSpecDetProc =movieData.getProcessIndex('SpeckleDetectionProcess',1,1);     
if isempty(iSpecDetProc)
    error(['Speckle detection has not yet been performed'...
    'on this movie! Please run first!!']);
end        
%Check that there is a valid output
specDetProc = movieData.processes_{iSpecDetProc};
% if ~specDetProc.checkChannelOutput(p.ChannelIndex)
%     error(['Each channel must have speckles! ' ...
%         'Please apply speckle detection to all needed channels before '...
%         'running flow analysis!'])
% end
p.MaskChannelIndex = specDetProc.funParams_.MaskChannelIndex;

% Create mask directory if several masks need to be merged
if length(p.MaskChannelIndex) >1
    %Get the indices of any previous mask intersection process
    iMaskProc = movieData.getProcessIndex('MaskIntersectionProcess',1,0);
else  
    iMaskProc =movieData.getProcessIndex('MaskRefinementProcess',1,1);
end

if isempty(iMaskProc)
    error(['Mask refinement has not yet been performed '...
        'on this movie! Please run first!!']);
end
%Check that there is a valid output
maskProc = movieData.processes_{iMaskProc};
% if ~maskProc.checkChannelOutput(p.ChannelIndex)
%     error(['Each channel must have masks !' ...
%         'Please apply mask refinement to all needed channels before'...
%         'running speckle tracking!'])
% end

% Test the presence and output validity of the speckle tracking process
% iFlowProc =movieData.getProcessIndex(p.FlowProcess,1,1);     
% if isempty(iFlowProc)
%     error([eval([p.FlowProcess '.getName']) ' has not yet been performed'...
%     ' on this movie! Please run first!!']);
% end        

%Check that there is a valid output
% flowProc = movieData.processes_{iFlowProc};
flowProc = movieData.processes_{iFlowProc};
% if ~flowProc.checkChannelOutput(p.ChannelIndex)
%     error(['Each channel must have flow! Please apply '...
%         eval([p.FlowProcess '.getName']) ' to all needed channels before '...
%         'running flow analysis!'])
% end
% this input part should be updated once the algorithm is completed.

%% Import the flow data per each frame - the flow should've been calculated
% between only the two adjacent frames
%% --------------- Kinetic analysi ---------------%%% 

disp('Starting reconstructing forces...')
% Format string for zero-padding file names
% Set up the input directories
nChan=numel(movieData.channels_);
outputDir=cell(1,nChan);
fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
numStr = @(frame) num2str(frame,fString);
% Anonymous functions for reading input/output
logMsg = @(chan) ['Please wait, analyzing flow for channel ' num2str(chan)];
outFile=@(chan,frame) [outputDir{chan} filesep 'intraForces_' numStr(frame) '.mat'];

speedMapLimits=cell(1,nChan);
flowLimits=cell(1,nChan);
channelLog=cell(1,numel(p.ChannelIndex));

for i=1:numel(p.ChannelIndex)
    iChan = p.ChannelIndex(i);
    % Log display
%     channelLog{i} = sprintf('Channel %g: %s\n',iChan,inFilePaths{1,iChan});
    disp(logMsg(iChan))
    
    if ishandle(wtBar), waitbar(0,wtBar,'Loading masks and images...'); end
%     maskNames = maskProc.getOutMaskFileNames(iChan);
%     inMask=@(frame) [flowAnProc.inFilePaths_{2,iChan}...
%         filesep maskNames{1}{frame}];
    mask=true([movieData.imSize_ nFrames]);
    stack=zeros([movieData.imSize_ nFrames]);
    for jj = 1:nFrames
%         mask(:,:,jj) = logical(imread(inMask(jj)));
        mask(:,:,jj) = maskProc.loadChannelOutput(iChan,jj);
        stack(:,:,jj) = movieData.channels_(iChan).loadImage(jj);
    end

    % Load candidates and generate Nx3 matrices with position and intensity
    % Replace fsmTrackFillSpeckleList
    if ishandle(wtBar), waitbar(0,wtBar,['Loading flow for channel ' num2str(iChan)']); end
    
    flow = flowProc.loadChannelOutput(iChan,'output','flow');
    flow=flow(1:end-1);
    % Remove NaNs 
    for jj=find(~cellfun(@isempty,flow))
        flow{jj}=flow{jj}(~isnan(flow{jj}(:,3)),:);
    end           

    % Interpolate field
    if ishandle(wtBar), waitbar(.25,wtBar,['Interpolating flow for channel ' num2str(iChan)']); end
    [Md,Ms,E,S,stats] =  analyzeFlow(flow,p.timeWindow+1,p.maxCorLength,...
        'noise',1,'error',1);
    
   %% Geometry detemination for meshing (decsg and initmesh)
    % It might be good to use the mask from mask process
    % and determine free boundary vs. inner boundary
    % by looking at its straightness
    % Take the binary image into geometry with 4 edges (3 inner + 1 free)
    % mask is mask(:,:,j)
    for jj = 1:nFrames
        % Adjust mask boundaries with speckle boundaries
        curFlow = flow{jj};
        LeftUpperCorner(1:2) = [min(curFlow(:,2)), min(curFlow(:,1))];
        RightLowerCorner(1:2) = [max(curFlow(:,2)), max(curFlow(:,1))];
        
        mask(RightLowerCorner(2)+1:end,:,jj) = 0; % for bottom
        mask(:,RightLowerCorner(1)+1:end,jj) = 0; % for right
        mask(1:LeftUpperCorner(2)-1,:,jj) = 0; % for top
        mask(:,1:LeftUpperCorner(1)-1,jj) = 0; % for left
 
%         B{1}(ymaxIdx,1) = RightLowerCorner(2); % for bottom
%         B{1}(xmaxIdx,2) = RightLowerCorner(1); % for right
%         B{1}(yminIdx,1) = LeftUpperCorner(2); % for top
%         B{1}(xminIdx,2) = LeftUpperCorner(1); % for left
        %% geometry
        % get the image boundaries
        B = bwboundaries(mask(:,:,jj));
        % Build a Geometry Description Matrix - I'll use curve (polygon
        % solid)
        % find the straight lines first
        xmin = min(B{1}(:,2));
        xmax = max(B{1}(:,2));
        ymin = min(B{1}(:,1));
        ymax = max(B{1}(:,1));
        minImgSize = 5; % edge length should be more than 5 pixel.
        nFreeEdge = 0; % the id number of free edge (1: left, 2:top, 3:right, 4:bottom)
        
        % rotating from left, top, right to bottom, check the straightness
        % (if the line is inner boundary vs. free boundary (curved))
        % for left edge
        xminIdx = B{1}(:,2)==xmin; % logical index
        if sum(xminIdx) > max(minImgSize, 0.1*(ymax-ymin)) && ...
            (sum(diff(sort(B{1}(xminIdx,1)))==1) == sum(xminIdx)-1  || ...
            sum(diff(sort(B{1}(xminIdx,1)))==1) == sum(xminIdx)-2) % see if they are consecutive
            curveL = [xmin max(B{1}(xminIdx,1)) xmin min(B{1}(xminIdx,1))]; % from bottom to top
            [~,curveLIdx1] = max(B{1}(xminIdx,1));
            [~,curveLIdx2] = min(B{1}(xminIdx,1));
            curveLIdxIdx = find(xminIdx);
            curveLIdx = curveLIdxIdx([curveLIdx1;curveLIdx2]);
        else
            curveL = [];
            xminIdx = false(size(xminIdx));
            nFreeEdge = 1;
        end
        % for top edge
        yminIdx = B{1}(:,1)==ymin; % logical index
        if sum(yminIdx) > max(minImgSize, 0.1*(xmax-xmin)) && ...
            (sum(diff(sort(B{1}(yminIdx,2)))==1) == sum(yminIdx)-1  || ...
            sum(diff(sort(B{1}(yminIdx,2)))==1) == sum(yminIdx)-2) % see if they are consecutive
            curveT = [min(B{1}(yminIdx,2)) ymin max(B{1}(yminIdx,2)) ymin]; % from left to right
            [~,curveTIdx1] = min(B{1}(yminIdx,2));
            [~,curveTIdx2] = max(B{1}(yminIdx,2));
            curveTIdxIdx = find(yminIdx);
            curveTIdx = curveTIdxIdx([curveTIdx1;curveTIdx2]);
        else
            % curve approximation with ...
            curveT = [];
            yminIdx = false(size(yminIdx));
            nFreeEdge = 2;
        end
        % for right edge
        xmaxIdx = B{1}(:,2)==xmax; % logical index
        if sum(xmaxIdx) > max(minImgSize, 0.1*(ymax-ymin)) && ...
            (sum(diff(sort(B{1}(xmaxIdx,1)))==1) == sum(xmaxIdx)-1  || ...
            sum(diff(sort(B{1}(xmaxIdx,1)))==1) == sum(xmaxIdx)-2) % see if they are consecutive
            curveR = [xmax min(B{1}(xmaxIdx,1)) xmax max(B{1}(xmaxIdx,1))]; % from top to bottom
            [~,curveRIdx1] = min(B{1}(xmaxIdx,1));
            [~,curveRIdx2] = max(B{1}(xmaxIdx,1));
            curveRIdxIdx = find(xmaxIdx);
            curveRIdx = curveRIdxIdx([curveRIdx1;curveRIdx2]);
        else
            curveR = [];
            xmaxIdx = false(size(xmaxIdx));
            nFreeEdge = 3;
        end
        % for bottom edge
        ymaxIdx = B{1}(:,1)==ymax; % logical index
        if sum(ymaxIdx) > max(minImgSize, 0.1*(xmax-xmin)) && ...
            (sum(diff(sort(B{1}(ymaxIdx,2)))==1) == sum(ymaxIdx)-1  || ...
            sum(diff(sort(B{1}(ymaxIdx,2)))==1) == sum(ymaxIdx)-2) % see if they are consecutive
            curveB = [max(B{1}(ymaxIdx,2)) ymax min(B{1}(ymaxIdx,2)) ymax]; % from bottom to top
            [~,curveBIdx1] = max(B{1}(ymaxIdx,2));
            [~,curveBIdx2] = min(B{1}(ymaxIdx,2));
            curveBIdxIdx = find(ymaxIdx);
            curveBIdx = curveBIdxIdx([curveBIdx1;curveBIdx2]);
        else
            curveB = [];
            ymaxIdx = false(size(ymaxIdx));
            nFreeEdge = 4;
        end
        % curve approximation with ...
        freeIdx = ~(xminIdx | yminIdx | xmaxIdx | ymaxIdx);% index for free edge
        freeIdx(max(find(freeIdx,1)-1,1)) = true;
        freeIdx(min(find(freeIdx,1,'last')+1,end)) = true;
        
        % determine distance between nodes based on speckle density
        numSpeckles=length(flow);
        areaImg=prod(RightLowerCorner-LeftUpperCorner);
        interSpecDist=ceil(sqrt(areaImg/numSpeckles));  % the number of skipping points for equi-spatial sampling of curves

        freeIdxIdx = find(freeIdx,1):interSpecDist:find(freeIdx,1,'last');
        if freeIdxIdx(end) ~= find(freeIdx,1,'last')
            freeIdxIdx(end) = find(freeIdx,1,'last');
        end
        
%         freeIdxS = false(size(freeIdx));
%         freeIdxS(freeIdxIdx) = true;
        
        switch nFreeEdge
            case 1
                curveLIdx = freeIdxIdx;
             case 2
                curveTIdx = freeIdxIdx;
            case 3
                curveRIdx = freeIdxIdx;
            case 4
                curveBIdx = freeIdxIdx;
            case 0
                disp('Something is wrong. There is no free edge. Check your boundary condition.')
        end
        
%         allEdgeIdx = xminIdx | yminIdx | xmaxIdx | ymaxIdx | freeIdxS;
        np = length(curveLIdx)+length(curveTIdx)+length(curveRIdx)+length(curveBIdx)-4; % the number of polygon segments
        gd = zeros(2+2*np,1); % initialization of geometry description matrix
        gd(1) = 2; % represents polygon solid
        gd(2) = np;
        gd(3:2+np) = B{1}([curveLIdx(1:end-1); curveTIdx(1:end-1)'; curveRIdx(1:end-1); curveBIdx(1:end-1)],2); % x-coordinates
        gd(3+np:2+2*np) = B{1}([curveLIdx(1:end-1); curveTIdx(1:end-1)'; curveRIdx(1:end-1); curveBIdx(1:end-1)],1); % y-coordinates
        dl = decsg(gd); % decompose constructive solid geometry into minimal regions
        %% mesh
        [p,e,t]=initmesh(dl,'hmax',2*interSpecDist); 
        curActin = movieData.getChannel(1).loadImage(jj);
        figure, imshow(curActin,[])
        hold on
        plot(curFlow(:,2),curFlow(:,1),'y.')
        quiver(curFlow(:,2),curFlow(:,1),3*(curFlow(:,4)-curFlow(:,2)),3*(curFlow(:,3)-curFlow(:,1)),0,'y')
        pdemesh(p,e,t)
        %% PDE definition for continuum mechanics, plane stress
        % Now I have a mesh information, with which I can solve for forward
        % solution for each basis function
        % I am using PDE structural mechanics, plane stress
        % governing equation:
        %           -div(c*grad(u)) = 0
        % where c11 = [2*G+mu 0; 0 G], c12 = [0 G; mu 0], c21=[0 G;mu 0],
        % c22 = [G 0; 0 2*G+mu];
        % Material properties
        E = 1e1; % Elastic modulus, Pa = N/m^2. I need to think about how to make it spatially heterogeneous later
        Nu = 0.29; % Poisson's ratio for actin network. need to find out more correct one
        G = E/(2*(1+Nu)); % Shear modulus, Pa
        mu = 2*G*Nu/(1-Nu); % mu
        c11 = [2*G+mu 0; 0 G];
        c12 = [0 G; mu 0];
        c21=[0 G;mu 0];
        c22 = [G 0; 0 2*G+mu];
        c2d = [c11 c12;c21 c22]; % 2Nx2N flattened coefficient c
        c = c2d(:); % vectorized c
        a = 0; %a=0 for structural mechanics
        f = [0 0]';
        %% Boundary condition definition
        nPDE = 2; % two dependent variables, ux and uy
        pb = pde(nPDE);
        % Create a geometry entity
        pg = pdeGeometryFromEdges(dl);
        % Set Diri
%         % Using speckle postions for mesh generation
%         dt=delaunayTriangulation(flow(:,2),flow(:,1));
%         triplot(dt)
        
        % Boundary condition (Free edge and Inner boundary edge)
        % Define the boundary condition vector, b, 
        % for the boundary condition u=x^2-y^2.
        % For each boundary segment, the boundary 
        
    end
    
    
    
    % 
    
end
    




























