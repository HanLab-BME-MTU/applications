function [M]=calcFwdMapFastBEM(x_vec_u, y_vec_u, forceMesh, E, meshPtsFwdSol, doPlot)

if nargin < 5
    meshPtsFwdSol=[];
end

if nargin < 6 || isempty(doPlot)
    doPlot=0;
end

% try to load the lookup table:
try
    load('basisClassTbl.mat');
catch
    basisClassTbl=struct([]) ;
end
addAtLeastOneToTbl=0;

forceSpan=1;
imgRows=1024;
imgCols=1344;
    

% transform to column vectors:
x_vec_u=x_vec_u(:);
y_vec_u=y_vec_u(:);


% test if displacement vectors have been measured only at integer
% positions:
diff_x_u=abs(x_vec_u-round(x_vec_u));
diff_y_u=abs(y_vec_u-round(y_vec_u));

% test if the basis function for the force are located only at integer
% positions:
allNodes=vertcat(forceMesh.basis(:).node);
diff_xy_f=abs(allNodes-round(allNodes));

if sum(diff_x_u(:))+sum(diff_y_u(:))<10^(-3) && sum(diff_xy_f(:))<10^(-3)
    method='direct';
else
    method='*cubic';
end

% determine the size of M
numBasis=length(forceMesh.basis);
numPts  =length(x_vec_u);

M=NaN*zeros(2*numPts,2*numBasis);

% here actually each basis function hast to be evaluated twice, which means
% it actually give two values:
ux=zeros(length(x_vec_u),2*numBasis);
uy=zeros(length(y_vec_u),2*numBasis);

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
    
% calculate the basis solution for all basis classes:
numClass=length(forceMesh.basisClass);
for class=1:numClass
    display(['Work on class: ',num2str(class),' of: ',num2str(numClass),'. Each class [~2x2min]:... ']);
    
    % Integration bounds used in the refine step, in case the fwdSolution
    % has to be calculated:
    xbd_min=min(forceMesh.basisClass(class).neighPos(:,1));
    xbd_max=max(forceMesh.basisClass(class).neighPos(:,1));
    ybd_min=min(forceMesh.basisClass(class).neighPos(:,2));
    ybd_max=max(forceMesh.basisClass(class).neighPos(:,2));

    % try to find the basis class in the table base:
    basisClassIn=forceMesh.basisClass(class);
    [idBestMatch]=findBasisClassInTbl(basisClassTbl,basisClassIn,xrangeReq,yrangeReq,meshPtsFwdSol);
    
    % now run through the x/y-comp. and either pull the fwdSol from the
    % table base or calculate it from scratch. The latter will be saved in
    % the basisClassTbl for later use
    for oneORtwo=1:2        
        if ~isempty(idBestMatch)
            % Then we can use the stored solution.
            % scale the basis solution with the right Youngs modulus. This
            % works for the boussinesq-BCs but might fail for more general BCs:
            scaleE = basisClassTbl(end).uSol.E/E; % for details see Landau Lifschitz p32.
            ux_model_pix = scaleE*double(basisClassTbl(idBestMatch).uSol.comp(oneORtwo).ux);
            uy_model_pix = scaleE*double(basisClassTbl(idBestMatch).uSol.comp(oneORtwo).uy);
            x_model_pix  = double(basisClassTbl(idBestMatch).uSol.x);
            y_model_pix  = double(basisClassTbl(idBestMatch).uSol.y);
        else
            display('Could not find a good match in the tablebase! Have to calculate the solution!')
            % no good match has been found in the table, we have to calculate the
            % solution from scratch. Here we could force xrangeSol>xrangeReq
            % and meshPtsFwdSol>meshPtsFwdSol_min, such that it is more
            % likely that this solution can be used in the future. We will
            % store this solution only for method='direct'!
            if forceSpan==1
                dxSol=max(imgCols,dx);
                dySol=max(imgRows,dy);
                
                xrangeSol=[-dxSol dxSol]';
                yrangeSol=[-dySol dySol]';
            else
                xrangeSol=xrangeReq;
                yrangeSol=yrangeReq;
            end               
            % calculate the solution:
            [ux_model, uy_model, x_model, y_model]=fwdSolution(xrangeSol,yrangeSol,E,xbd_min,xbd_max,ybd_min,ybd_max,forceMesh.basisClass(class).basisFunc(oneORtwo).f_intp_x,forceMesh.basisClass(class).basisFunc(oneORtwo).f_intp_y,'fft','noIntp',meshPtsFwdSol);
            % check if the sampling is fine enough for method 'direct':
            if strcmp(method,'direct')
                % This works perfectly for all mesh types as long as the
                % displacment and force nodes are defined at integer positions!
                x_spacing=x_model(2,2)-x_model(1,1);
                y_spacing=y_model(2,2)-y_model(1,1);
                if x_spacing<=1 && y_spacing<=1
                    % Only if the spacing is <1 we have oversampled, interpolate to
                    % integer positions:
                    x_model_pix=x_model(1,1):1:x_model(end,end);
                    y_model_pix=y_model(1,1):1:y_model(end,end);
                    
                    [x_model_pix y_model_pix]=meshgrid(x_model_pix,y_model_pix);
                    
                    %interpolate the solution to the integer positions:
                    ux_model_pix= interp2(x_model, y_model, ux_model, x_model_pix, y_model_pix, '*cubic');  %This is ux(:,j)
                    uy_model_pix= interp2(x_model, y_model, uy_model, x_model_pix, y_model_pix, '*cubic');  %This is uy(:,j)
                else
                    display('Have switched over to *cubic. But is this really necessary? It might well be that even if undersampled the upper search will produce the same result as an interpolation')
                    method='*cubic';
                end
            end
        end
        
        toDoBasis=find(vertcat(forceMesh.basis.class)==class)';
        lgthToDoBasis=length(toDoBasis);
        display(['Evaluate ',num2str(lgthToDoBasis),' basis functions'])
        for basisID=toDoBasis
            % lgthToDoBasis=length(toDoBasis);
            % displayText=[num2str(basisID),' of ',num2str(lgthToDoBasis)];
            % progressText(basisID/lgthToDoBasis,displayText);
            % Interpolate the basis-solution:
            xShift = forceMesh.basis(basisID).node(1);
            yShift = forceMesh.basis(basisID).node(2);
            
            if strcmp(method,'direct')
                % Instead of interpolation we can simply search for the
                % correct values in the matrix:
                xmin_pix=x_model_pix(1,1)+xShift;
                ymin_pix=y_model_pix(1,1)+yShift;
                
                % shift the x,y values to indices:
                px=x_vec_u-xmin_pix+1; % +1 because the matrix index starts at 1.
                py=y_vec_u-ymin_pix+1; % +1 because the matrix index starts at 1.
                
                % get the indices:
                [ptInd]=sub2ind(size(x_model_pix),py,px);
                
                if oneORtwo==1
                    M(1:numPts    ,basisID)          = ux_model_pix(ptInd);
                    M(numPts+1:end,basisID)          = uy_model_pix(ptInd);
                elseif oneORtwo==2
                    M(1:numPts    ,basisID+numBasis) = ux_model_pix(ptInd);
                    M(numPts+1:end,basisID+numBasis) = uy_model_pix(ptInd);
                end
                
%                 % This might be a bit more robust but is slower:
%                 forceNearest=0;
%                 if forceNearest==1
%                     if oneORtwo==1
%                         % Then the interpolants of the first function are:
%                         M(1:numPts    ,basisID)          = interp2(x_model_pix+xShift, y_model_pix+yShift, ux_model_pix, x_vec_u, y_vec_u, '*nearest');  %This is ux(:,j)
%                         M(numPts+1:end,basisID)          = interp2(x_model_pix+xShift, y_model_pix+yShift, uy_model_pix, x_vec_u, y_vec_u, '*nearest');  %This is uy(:,j)
%                     elseif oneORtwo==2
%                         % Then the interpolants of the second function are:  (:,j+numBasis)
%                         M(1:numPts    ,basisID+numBasis) = interp2(x_model_pix+xShift, y_model_pix+yShift, ux_model_pix, x_vec_u, y_vec_u, '*nearest');  %This is ux(:,j+numBasis)
%                         M(numPts+1:end,basisID+numBasis) = interp2(x_model_pix+xShift, y_model_pix+yShift, uy_model_pix, x_vec_u, y_vec_u, '*nearest');  %This is uy(:,j+numBasis)
%                     end
%                 end
            else
                if oneORtwo==1
                    % Then the interpolants of the first function are:
                    M(1:numPts    ,basisID)          = interp2(x_model+xShift, y_model+yShift, ux_model, x_vec_u, y_vec_u, method);  %This is ux(:,j)
                    M(numPts+1:end,basisID)          = interp2(x_model+xShift, y_model+yShift, uy_model, x_vec_u, y_vec_u, method);  %This is uy(:,j)
                elseif oneORtwo==2
                    % Then the interpolants of the second function are:  (:,j+numBasis)
                    M(1:numPts    ,basisID+numBasis) = interp2(x_model+xShift, y_model+yShift, ux_model, x_vec_u, y_vec_u, method);  %This is ux(:,j+numBasis)
                    M(numPts+1:end,basisID+numBasis) = interp2(x_model+xShift, y_model+yShift, uy_model, x_vec_u, y_vec_u, method);  %This is uy(:,j+numBasis)
                end
            end
        end
        if isempty(idBestMatch) && strcmp(method,'direct')
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
            basisClassTbl(end).uSol.comp(oneORtwo).ux = single(ux_model_pix*E); % the factor E is to scale the solution to u(E=1,f)
            basisClassTbl(end).uSol.comp(oneORtwo).uy = single(uy_model_pix*E); % the factor E is to scale the solution to u(E=1,f)
            basisClassTbl(end).uSol.x = int16(x_model_pix);
            basisClassTbl(end).uSol.y = int16(y_model_pix);
            
            % enter parameters:
            basisClassTbl(end).uSol.xrange       = xrangeSol;
            basisClassTbl(end).uSol.yrange       = yrangeSol;
            basisClassTbl(end).uSol.E            = 1; % this could be more general!
            basisClassTbl(end).uSol.method       ='fft';
            basisClassTbl(end).uSol.meshPtsFwdSol= meshPtsFwdSol;
            addAtLeastOneToTbl=1;
        end
    end    
end
% save the new table base
if addAtLeastOneToTbl
    save('basisClassTbl.mat', 'basisClassTbl', '-v7.3');
end

% plot an example to see if it works correctly
if doPlot==1
    ind=1;
    if forceMesh.numBasis>ind-1
        xmin=min(x_vec_u);
        ymin=min(y_vec_u);
        xmax=max(x_vec_u);
        ymax=max(y_vec_u);
        figure(11)
        quiver(x_vec_u,y_vec_u,ux(:,ind),uy(:,ind))
        hold on
        quiver(x_vec_u,y_vec_u,ux(:,ind+forceMesh.numBasis),uy(:,ind+forceMesh.numBasis))
        xlim([xmin xmax])
        ylim([ymin ymax])
        hold off
    end
end