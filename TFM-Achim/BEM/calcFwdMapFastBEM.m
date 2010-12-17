function [M]=calcFwdMapFastBEM(x_vec_u, y_vec_u, forceMesh, E, meshPtsFwdSol, doPlot)

if nargin < 5
    meshPtsFwdSol=[];
end

if nargin < 6 || isempty(doPlot)
    doPlot=0;
end

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

% calculate the basis solution for all basis classes:
numClass=length(forceMesh.basisClass);
for class=1:numClass
    display(['Work on class: ',num2str(class),' of: ',num2str(numClass),'. Each class [~2x2min]:... ']);
    % To make sure that the range over which the solution is calculated,
    % take the double of the initial x and y ranges:
    xmin=min(x_vec_u); xmax=max(x_vec_u);
    ymin=min(y_vec_u); ymax=max(y_vec_u);
    dx=xmax-xmin;
    dy=ymax-ymin;
    
    xrange=[-dx dx]';
    yrange=[-dy dy]';
    
    % Integration bounds used in the refine step:
    xbd_min=min(forceMesh.basisClass(class).neighPos(:,1));
    xbd_max=max(forceMesh.basisClass(class).neighPos(:,1));
    ybd_min=min(forceMesh.basisClass(class).neighPos(:,2));
    ybd_max=max(forceMesh.basisClass(class).neighPos(:,2));
    
    for oneORtwo=1:2
        % calculate basis solution:
        [ux_model, uy_model, x_model, y_model]=fwdSolution(xrange,yrange,E,xbd_min,xbd_max,ybd_min,ybd_max,forceMesh.basisClass(class).basisFunc(oneORtwo).f_intp_x,forceMesh.basisClass(class).basisFunc(oneORtwo).f_intp_y,'fft','noIntp',meshPtsFwdSol);
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
        
        
        for basisID=1:numBasis
            displayText=[num2str(basisID),' of ',num2str(numBasis)];
            progressText(basisID/numBasis,displayText);
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
    end
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