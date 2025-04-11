function [M]=calcFwdMapFEM(x_vec_u, y_vec_u, forceMesh, E,varargin)
% Synoposis:  M=calcFwdMapFEM(x_vec_u, y_vec_u, forceMesh, E)

ip = inputParser;
ip.addRequired('x_vec_u',@isnumeric);
ip.addRequired('y_vec_u',@isnumeric);
ip.addRequired('forceMesh',@isstruct);
ip.addRequired('E',@isscalar);
ip.addOptional('meshPtsFwdSol',[],@isscalar);
ip.addOptional('doPlot',0,@isscalar);
ip.addParameter('basisClassTblPath','basisClassTbl.mat',@ischar);
ip.addParameter('wtBar',-1,@isscalar);
ip.addParameter('imgRows',1024,@isscalar);
ip.addParameter('imgCols',1334,@isscalar);
ip.addParameter('thickness',472,@isscalar); % default assuming 34 um with 72 nm/pix resolution
ip.addParameter('PoissonRatio',0.5,@isscalar); 
ip.parse(x_vec_u, y_vec_u, forceMesh, E,varargin{:})
meshPtsFwdSol=ip.Results.meshPtsFwdSol;
doPlot=ip.Results.doPlot;
basisClassTblPath=ip.Results.basisClassTblPath;
wtBar=ip.Results.wtBar;
imgRows = ip.Results.imgRows;
imgCols = ip.Results.imgCols;
thickness = ip.Results.thickness;    
v = ip.Results.PoissonRatio;

disp('calcFwdMapFEM')

% try to load the lookup table:
try
    basisClassTblData=load(basisClassTblPath);
    basisClassTbl=basisClassTblData.basisClassTbl;
    addAtLeastOneToTbl=0;
catch %ME
    basisClassTbl=struct([]) ;
    addAtLeastOneToTbl=1;
end


% for a field of view with 1024x1344 pix and a square grid with grid size
% of 10 pix, the average rel. difference between the solutions with 2^11 and
% 2^12 pts for the fwd solution is less than 0.27% for all positions with
% non-vanishing forces. The maximum error is 4%. These results indicate
% that actually 2^11 points are fine for calculating the basis solutions
% for such a configuration. 2^12 should be fine for sure!
% forceSpan=1;
% imgRows=1024;
% imgCols=1344;
    

% transform to column vectors:
x_vec_u=x_vec_u(:);
y_vec_u=y_vec_u(:);


% test if displacement vectors have been measured only at integer
% positions:
diff_x_u = abs(x_vec_u-round(x_vec_u));
diff_y_u = abs(y_vec_u-round(y_vec_u));

% test if the basis function for the force are located only at integer
% positions:
allNodes = vertcat(forceMesh.basis(:).node);
diff_xy_f = abs(allNodes-round(allNodes));

if sum(diff_x_u(:))+sum(diff_y_u(:))<10^(-3) && sum(diff_xy_f(:))<10^(-3)
    method='direct';
else
    method='*cubic';
end

method = '*cubic';

% determine the size of M
numBasis=length(forceMesh.basis);
numPts  =length(x_vec_u);

% Initialize M
M=NaN(2*numPts,2*numBasis);

% here actually each basis function hast to be evaluated twice, which means
% it actually give two values:
%ux=zeros(length(x_vec_u),2*numBasis);
%uy=zeros(length(y_vec_u),2*numBasis);

% To make sure that the range over which the solution is calculated,
% take the double of the initial x and y ranges:
xmin=min(x_vec_u); xmax=max(x_vec_u);
ymin=min(y_vec_u); ymax=max(y_vec_u);
dx=xmax-xmin;
dy=ymax-ymin;

% The minimum x/y-range over which the basis solution has/had to be
% calculated:
xrangeReq=[-dx dx]';
yrangeReq=[-dy dy]';

dxSol=max(imgCols,dx);
dySol=max(imgRows,dy);

xrangeSol=[-dxSol dxSol]';
yrangeSol=[-dySol dySol]';

%[nodePtsX, nodePtsY] = meshgrid(linspace(min(allNodes(:,1)),max(allNodes(:,1)),numel(unique(allNodes(:,1)))),linspace(min(allNodes(:,1)),max(allNodes(:,1)),numel(unique(allNodes(:,1)))));

grooveWidth = 5; %should be 1/2 of the groove width that is given in fwdSolution OR MAYBE NOT??
fprintf('calcFwdMapFEM groove width: %1.1f micron \n',grooveWidth/10);
%numPix_x = abs(xrangeSol(1))+xrangeSol(2);
halfSide = xrangeSol(2);
%numGrooves = floor(ceil(numPix_x / grooveWidth)/2); %determine number of grooves
% edgePad = mod(numPix_x,grooveWidth)/2; %calculate dead space at edges of substrate where no full groove can fit
% grooveEdges = -halfSide+edgePad-grooveWidth/2:grooveWidth:halfSide-edgePad+grooveWidth/2; %calculate groove edges in x-direction
% if ~mod(numel(grooveEdges), 2) == 0 && grooveWidth ~= 8 %if number of groove edges is odd we have to refine
%     grooveEdges = -halfSide+edgePad:grooveWidth:halfSide-edgePad;
%     numGrooves = numGrooves - 1;
% end
% if grooveWidth == 20
%     grooveEdges = grooveEdges(3:end-2);
%     numGrooves = numGrooves - 1;
% elseif grooveWidth == 8
%     grooveEdges = grooveEdges(3:end-2);
%     numGrooves = numGrooves - 1;
% else
%     grooveEdges = grooveEdges(2:end-1); %ensure a groove is placed directly in the center of the substrate
% end
% 
% %calculate number of grooves touched by force area
% numPreservedGrooves = 1; %desired number of extra grooves adjacent to outermost force-touched groove
% if grooveWidth / 2 <= xmax || grooveWidth / 2 <= ymax
%     overlap = ceil((xmax / grooveWidth)) / 2 - 1;
%     numPreservedGrooves = overlap + numPreservedGrooves;
% end
% %numPreservedGrooves must always be 2 or a multiple of 2
% numPreservedGrooves = numPreservedGrooves * 2; 
% if numPreservedGrooves == 1
%     numPreservedGrooves = numPreservedGrooves + 1;
% end
% %trim outer grooves beyond substrate area
% if ~grooveEdges(1) < xrangeSol(1)
%     grooveEdges = grooveEdges(length(grooveEdges)/2 - numPreservedGrooves:length(grooveEdges)/2 + 1 + numPreservedGrooves);
%     numGrooves = length(grooveEdges) / 2;
% end

grooveEdges = [-grooveWidth/2, grooveWidth/2];
grooveFront = grooveWidth/2;
numGrooves = 1;
while grooveFront < dx/2
    grooveEdges = [-grooveFront - 2*grooveWidth, -grooveFront - grooveWidth, grooveEdges, grooveFront + grooveWidth, grooveFront + 2*grooveWidth]; %#ok<*AGROW>
    numGrooves = numGrooves + 2;
    grooveFront = grooveFront + 2*grooveWidth;
end
grooveEdges = grooveEdges + halfSide/2;

% grooveSet = flip(grooveEdges(1:length(grooveEdges)/2));
% grooveSet = flip([grooveSet; grooveEdges(length(grooveEdges)/2+1:end)]);

ptsOnGrooves = false(size(allNodes(:,1)));
groovePolysUnfuzzed{numGrooves} = []; g = 1;
for i = 1:2:length(grooveEdges)
    groovePolysUnfuzzed{g} = polyshape([grooveEdges(i), grooveEdges(i), grooveEdges(i+1), grooveEdges(i+1)],[halfSide, -halfSide, -halfSide, halfSide]);
    xv = groovePolysUnfuzzed{g}.Vertices(:,1);
    yv = groovePolysUnfuzzed{g}.Vertices(:,2);
    ptsOnGrooves = ptsOnGrooves | inpolygon(allNodes(:,1),allNodes(:,2),xv,yv);
    g = g + 1;
end
ptsInGrooves = ~ptsOnGrooves;

% temp = ptsInGrooves;
% ptsInGrooves = ptsOnGrooves;
% ptsOnGrooves = temp;

if doPlot
figure,
scatter(allNodes(ptsOnGrooves,1),allNodes(ptsOnGrooves,2),'ko')
hold on
scatter(allNodes(ptsInGrooves,1),allNodes(ptsInGrooves,2),'b.')
hold off
end

disp('Begin new multi-geometric fwdsol test')

class = 1;
xbd_min=min(forceMesh.basisClass(class).neighPos(:,1));
xbd_max=max(forceMesh.basisClass(class).neighPos(:,1));
ybd_min=min(forceMesh.basisClass(class).neighPos(:,2));
ybd_max=max(forceMesh.basisClass(class).neighPos(:,2));

% try to find the basis class in the table base:
basisClassIn=forceMesh.basisClass(class);
[idBestMatch]=findBasisClassInTbl(basisClassTbl,basisClassIn,xrangeReq,yrangeReq,meshPtsFwdSol);

for oneORtwo = 1:2
    if ~isempty(idBestMatch)
        % Then we can use the stored solution.
        % scale the basis solution with the right Youngs modulus. This
        % works for the boussinesq-BCs but might fail for more general BCs:
        scaleE = basisClassTbl(end).uSol.E/E; % for details see Landau Lifschitz p32.
        %GROOVE TOP INFO
        ux_model_pix = scaleE*double(basisClassTbl(idBestMatch).uSol.comp(oneORtwo).ux);
        uy_model_pix = scaleE*double(basisClassTbl(idBestMatch).uSol.comp(oneORtwo).uy);
        x_model_pix  = double(basisClassTbl(idBestMatch).uSol.x);
        y_model_pix  = double(basisClassTbl(idBestMatch).uSol.y);
        %GROOVE VALLEY INFO
        try
            ux_model_pix_flat = scaleE*double(basisClassTbl(idBestMatch).uSol.comp(oneORtwo).flatux);
            uy_model_pix_flat = scaleE*double(basisClassTbl(idBestMatch).uSol.comp(oneORtwo).flatuy);
            x_model_pix_flat  = double(basisClassTbl(idBestMatch).uSol.flatx);
            y_model_pix_flat  = double(basisClassTbl(idBestMatch).uSol.flaty);
        catch ME
            msg = 'Saved basis class from non-FEM solution, rerun an FEM solution to generate an appropriate basis class.';
            causeException = MException('MATLAB:myCode:basisclassMismatch',msg);
            ME = addCause(ME,causeException);
            rethrow(ME)
        end 
    else
        [ux_model, uy_model, x_model, y_model]=fwdSolution(xrangeSol,yrangeSol,E,...
            xbd_min,xbd_max,ybd_min,ybd_max,forceMesh.basisClass(class).basisFunc(oneORtwo).f_intp_x,forceMesh.basisClass(class).basisFunc(oneORtwo).f_intp_y,...
            'FEM','noIntp',meshPtsFwdSol,thickness,v,false,false);
        
        if strcmp(method,'direct') || strcmp(method,'*cubic')
            % This works perfectly for all mesh types as long as the
            % displacment and force nodes are defined at integer positions!
            x_spacing=x_model(2,2)-x_model(1,1);
            y_spacing=y_model(2,2)-y_model(1,1);
            if x_spacing<=1 && y_spacing<=1
                % Only if the spacing is <1 we have oversampled, interpolate to
                % integer positions:
                x_model_pix=x_model(1,1):1:x_model(end,end);
                y_model_pix=y_model(1,1):1:y_model(end,end);
                
                [x_model_pix,y_model_pix]=meshgrid(x_model_pix,y_model_pix);
                
                %interpolate the solution to the integer positions:
                ux_model_pix= interp2(x_model, y_model, ux_model, x_model_pix, y_model_pix); %, 'direct'); There is no such thing as direct, but only 'linear'  %This is ux(:,j)
                uy_model_pix= interp2(x_model, y_model, uy_model, x_model_pix, y_model_pix); %, 'direct');  %This is uy(:,j)
            else
                disp('Have switched over to *cubic. But is this really necessary? It might well be that even if undersampled the upper search will produce the same result as an interpolation')
                method='*cubic';
                pizInterval_x = round(x_spacing);
                pizInterval_y = round(y_spacing);
                x_model_pix=x_model(1,1):pizInterval_x:x_model(end,end);
                y_model_pix=y_model(1,1):pizInterval_y:y_model(end,end);
                
                [x_model_pix,y_model_pix]=meshgrid(x_model_pix,y_model_pix);
                
                %interpolate the solution to the integer positions:
                ux_model_pix= interp2(x_model, y_model, ux_model, x_model_pix, y_model_pix, 'cubic');  %This is ux(:,j)
                uy_model_pix= interp2(x_model, y_model, uy_model, x_model_pix, y_model_pix, 'cubic');  %This is uy(:,j)
            end
        end
        
        grooveHeight = 0;
        [ux_model_flat, uy_model_flat, x_model_flat, y_model_flat]=fwdSolutionFEM(grooveHeight,xrangeSol,yrangeSol,E,...
            xbd_min,xbd_max,ybd_min,ybd_max,forceMesh.basisClass(class).basisFunc(oneORtwo).f_intp_x,forceMesh.basisClass(class).basisFunc(oneORtwo).f_intp_y,...
            'FEM','noIntp',meshPtsFwdSol,thickness,v,false,false);
        
        if strcmp(method,'direct') || strcmp(method,'*cubic')
            % This works perfectly for all mesh types as long as the
            % displacment and force nodes are defined at integer positions!
            x_spacing_flat=x_model_flat(2,2)-x_model_flat(1,1);
            y_spacing_flat=y_model_flat(2,2)-y_model_flat(1,1);
            if x_spacing_flat<=1 && y_spacing_flat<=1
                % Only if the spacing is <1 we have oversampled, interpolate to
                % integer positions:
                x_model_pix_flat=x_model_flat(1,1):1:x_model_flat(end,end);
                y_model_pix_flat=y_model_flat(1,1):1:y_model_flat(end,end);
                
                [x_model_pix_flat,y_model_pix_flat]=meshgrid(x_model_pix_flat,y_model_pix_flat);
                
                %interpolate the solution to the integer positions:
                ux_model_pix_flat= interp2(x_model_flat, y_model_flat, ux_model_flat, x_model_pix_flat, y_model_pix_flat); %, 'direct'); There is no such thing as direct, but only 'linear'  %This is ux(:,j)
                uy_model_pix_flat= interp2(x_model_flat, y_model_flat, uy_model_flat, x_model_pix_flat, y_model_pix_flat); %, 'direct');  %This is uy(:,j)
            else
                disp('Have switched over to *cubic. But is this really necessary? It might well be that even if undersampled the upper search will produce the same result as an interpolation')
                method='*cubic';
                pizInterval_x_flat = round(x_spacing_flat);
                pizInterval_y_flat = round(y_spacing_flat);
                x_model_pix_flat=x_model_flat(1,1):pizInterval_x_flat:x_model_flat(end,end);
                y_model_pix_flat=y_model_flat(1,1):pizInterval_y_flat:y_model_flat(end,end);
                
                [x_model_pix_flat,y_model_pix_flat]=meshgrid(x_model_pix_flat,y_model_pix_flat);
                
                %interpolate the solution to the integer positions:
                ux_model_pix_flat= interp2(x_model_flat, y_model_flat, ux_model_flat, x_model_pix_flat, y_model_pix_flat, 'cubic');  %This is ux(:,j)
                uy_model_pix_flat= interp2(x_model_flat, y_model_flat, uy_model_flat, x_model_pix_flat, y_model_pix_flat, 'cubic');  %This is uy(:,j)
            end
        end
    end
    
    toDoBasis=find(vertcat(forceMesh.basis.class)==class)';
    lgthToDoBasis=length(toDoBasis);
    disp(['Evaluate ',num2str(lgthToDoBasis),' basis functions'])
    
    
    logMsg = 'Please wait, interpolating basis solutions';
    timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];
    tic;
    if ishandle(wtBar)
        wtBar = waitbar(0,wtBar,logMsg);
%         elseif feature('ShowFigureWindows'),
%             wtBar = waitbar(0,logMsg);
    end
    
    x_spacing=x_model_pix(2,2)-x_model_pix(1,1);
    y_spacing=y_model_pix(2,2)-y_model_pix(1,1);
    if x_spacing>1 || y_spacing>1
        method = '*cubic';
    end
    
    collatedxshifts = [];
    collatedyshifts = [];

    for i=1:numel(toDoBasis)
        basisID=toDoBasis(i);
        % lgthToDoBasis=length(toDoBasis);
        % displayText=[num2str(basisID),' of ',num2str(lgthToDoBasis)];
        % progressText(basisID/lgthToDoBasis,displayText);
        % Interpolate the basis-solution:
        xShift = forceMesh.basis(basisID).node(1);
        collatedxshifts = [collatedxshifts, xShift];
        yShift = forceMesh.basis(basisID).node(2);
        collatedyshifts = [collatedyshifts, yShift];
        if ptsInGrooves(i)
            if oneORtwo==1
                % Then the interpolants of the first function are:
                M(1:numPts    ,basisID)          = interp2(x_model_pix_flat+xShift, y_model_pix_flat+yShift, ux_model_pix_flat, x_vec_u, y_vec_u, method);  %This is ux(:,j)
                M(numPts+1:end,basisID)          = interp2(x_model_pix_flat+xShift, y_model_pix_flat+yShift, uy_model_pix_flat, x_vec_u, y_vec_u, method);  %This is uy(:,j)
            elseif oneORtwo==2
                % Then the interpolants of the second function are:  (:,j+numBasis)
                M(1:numPts    ,basisID+numBasis) = interp2(x_model_pix_flat+xShift, y_model_pix_flat+yShift, ux_model_pix_flat, x_vec_u, y_vec_u, method);  %This is ux(:,j+numBasis)
                M(numPts+1:end,basisID+numBasis) = interp2(x_model_pix_flat+xShift, y_model_pix_flat+yShift, uy_model_pix_flat, x_vec_u, y_vec_u, method);  %This is uy(:,j+numBasis)
            end
        elseif ptsOnGrooves(i)
            if oneORtwo==1
                % Then the interpolants of the first function are:
                M(1:numPts    ,basisID)          = interp2(x_model_pix+xShift, y_model_pix+yShift, ux_model_pix, x_vec_u, y_vec_u, method);  %This is ux(:,j)
                M(numPts+1:end,basisID)          = interp2(x_model_pix+xShift, y_model_pix+yShift, uy_model_pix, x_vec_u, y_vec_u, method);  %This is uy(:,j)
            elseif oneORtwo==2
                % Then the interpolants of the second function are:  (:,j+numBasis)
                M(1:numPts    ,basisID+numBasis) = interp2(x_model_pix+xShift, y_model_pix+yShift, ux_model_pix, x_vec_u, y_vec_u, method);  %This is ux(:,j+numBasis)
                M(numPts+1:end,basisID+numBasis) = interp2(x_model_pix+xShift, y_model_pix+yShift, uy_model_pix, x_vec_u, y_vec_u, method);  %This is uy(:,j+numBasis)
            end
        end
        
        % Update the waitbar
        if mod(i,5)==1 && ishandle(wtBar)
            ti=toc;
            waitbar(i/lgthToDoBasis,wtBar,...
                sprintf([logMsg timeMsg(ti*lgthToDoBasis/i-ti)]));
        end
        
    end

    if isempty(idBestMatch) %&& strcmp(method,'direct')
        % Then we have either calculated a previously unknown basis
        % Solution, or we have improved one (by increasing the range or
        % by increasing the meshPtsFwdSol). Enter a new entry only if
        % we are working on the x-comp. The y-comp will be treated in
        % the next loop and will be sorted in into the same
        % basisClassTbl-id:
        if oneORtwo==1
            numClassTbl=length(basisClassTbl);
            currBasisClass=forceMesh.basisClass(class);
            % strip of the basisFunc entry. We don't need this
            % information!
            basisClassTbl(numClassTbl+1).centerPos  = currBasisClass.centerPos;
            basisClassTbl(numClassTbl+1).numNeigh   = currBasisClass.numNeigh;
            basisClassTbl(numClassTbl+1).neighPos   = currBasisClass.neighPos;
            basisClassTbl(numClassTbl+1).dtBaseSup  = currBasisClass.dtBaseSup;
            basisClassTbl(numClassTbl+1).unitVolume = currBasisClass.unitVolume;
        end
        % enter the basis solutions:
        % Scale the basis solution with the right Youngs modulus. This
        % works for the boussinesq-BCs but might fail for more general BCs:
        % To store the forward solution, single precision should be
        % sufficient. x/y positions are integer anyways, store them in
        % int16 format:
        %GROOVE TOP BASIS INFO
        basisClassTbl(end).uSol.comp(oneORtwo).ux = single(ux_model_pix*E); % the factor E is to scale the solution to u(E=1,f)
        basisClassTbl(end).uSol.comp(oneORtwo).uy = single(uy_model_pix*E); % the factor E is to scale the solution to u(E=1,f)
        basisClassTbl(end).uSol.x = int16(x_model_pix);
        basisClassTbl(end).uSol.y = int16(y_model_pix);
        %GROOVE BOTTOM BASIS INFO
        basisClassTbl(end).uSol.comp(oneORtwo).flatux = single(ux_model_pix_flat*E); % the factor E is to scale the solution to u(E=1,f)
        basisClassTbl(end).uSol.comp(oneORtwo).flatuy = single(uy_model_pix_flat*E); % the factor E is to scale the solution to u(E=1,f)
        basisClassTbl(end).uSol.flatx = int16(x_model_pix_flat);
        basisClassTbl(end).uSol.flaty = int16(y_model_pix_flat);
        
        % enter parameters:
        basisClassTbl(end).uSol.xrange       = xrangeSol;
        basisClassTbl(end).uSol.yrange       = yrangeSol;
        basisClassTbl(end).uSol.E                = 1; % this could be more general!
        basisClassTbl(end).uSol.method       ='fft';
        basisClassTbl(end).uSol.meshPtsFwdSol= meshPtsFwdSol;
        basisClassTbl(end).uSol.gelHeight   = thickness; % this could be more general!
        addAtLeastOneToTbl=1;
    end
end

if addAtLeastOneToTbl
    save(basisClassTblPath, 'basisClassTbl','-v7.3');
end

% plot an example to see if it works correctly
if doPlot==1
    ind=1;
    if forceMesh.numBasis>ind-1
        xmin=min(x_vec_u);
        ymin=min(y_vec_u);
        xmax=max(x_vec_u);
        ymax=max(y_vec_u);
        ux = basisClassTbl(end).uSol.comp(oneORtwo).ux;
        uy = basisClassTbl(end).uSol.comp(oneORtwo).uy;
        figure(11)
        quiver(x_vec_u,y_vec_u,ux(:,ind),uy(:,ind))
        hold on
        quiver(x_vec_u,y_vec_u,ux(:,ind+forceMesh.numBasis),uy(:,ind+forceMesh.numBasis))
        xlim([xmin xmax])
        ylim([ymin ymax])
        hold off
    end
end