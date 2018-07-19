function [refinedForceField] = refineForceField(forceField,frameList,x_out,y_out)

if nargin < 2 || isempty(frameList)
    toDoList=(1:length(forceField));
else
    toDoList=frameList(:)';
end

if nargin < 3
    x_out=[];
    y_out=[];
end

ROI=[];
for frame=toDoList

    dPix=50;
    min_x=min(forceField(frame).pos(:,1))-dPix;
    max_x=max(forceField(frame).pos(:,1))+dPix;
    min_y=min(forceField(frame).pos(:,2))-dPix;
    max_y=max(forceField(frame).pos(:,2))+dPix;
    maxForcePlot=1/forceField(frame).par.gridSpacing*max(sqrt(forceField(frame).vec(:,1).^2+forceField(frame).vec(:,2).^2));

    figure(1)
    quiver(forceField(frame).pos(:,1),forceField(frame).pos(:,2),forceField(frame).vec(:,1)/maxForcePlot,forceField(frame).vec(:,2)/maxForcePlot,0,'r');
    title(['Force field no: ',num2str(frame)])
    axis equal
    xlim([min_x max_x])
    ylim([min_y max_y])
    set(gca,'YDir','reverse')

    if ~isempty(ROI)
        redraw=input('Do you want to re-draw the ROIs Y or [N]?: ','s');
    else
        redraw='yes';
    end
    
    if strcmp(redraw,'Y') || strcmp(redraw,'y') || strcmp(redraw,'yes') || strcmp(redraw,'Yes') || strcmp(redraw,'1')
        numROI=input('How many regions do you want to refine?: ');
        for j=1:numROI
            ROI{j}.curve=[];
            while size(ROI{j}.curve,1)<2
                display('Encircle the region of interest!')
                polygonObject = impoly(gca);
                ROI{j}.curve  = round(getPosition(polygonObject));
            end
            display('Thanks, the region has been acquired!...')
            % get all points of the forceMesh that are in the drawn regions:
            ptsX=forceField(1).par.forceMesh.p(:,1);
            ptsY=forceField(1).par.forceMesh.p(:,2);
            checkVec = inpolygon(ptsX,ptsY,ROI{j}.curve(:,1),ROI{j}.curve(:,2));    
            ROI{j}.ptsToRef=horzcat(ptsX(checkVec),ptsY(checkVec));
        end
        % These are all points from all regions:
        allPtsToRef=[];
        for j=1:numROI
            allPtsToRef=vertcat(allPtsToRef,ROI{j}.ptsToRef);
        end
        % sort out double entries:
        allPtsToRef=sortrows(allPtsToRef);
        allPtsToRef=removeDoublePoints(allPtsToRef);
    end

    % read out the important values:
    forceMesh    =forceField(frame).par.forceMesh;
    M_old        =forceField(frame).par.M;
    x            =forceField(frame).par.pos(:,1);
    y            =forceField(frame).par.pos(:,2);
    ux           =forceField(frame).par.u(1:end/2);
    uy           =forceField(frame).par.u(end/2+1:end);
    E            =forceField(frame).par.yModu_Pa;
    L            =forceField(frame).par.regParam;
    meshPtsFwdSol=forceField(frame).par.meshPtsFwdSol;



    doPlot=1;
    [refinedForceMesh]=refineMeshAndBasis(forceMesh,allPtsToRef,doPlot);

    %[fx, fy, x_out, y_out, M, pos_u, u, sol_coef]
    [fx, fy, x_out, y_out, M,     ~, ~, sol_coef] = refine_BEM_force_reconstruction(x,y,ux,uy,M_old,refinedForceMesh,E,L,meshPtsFwdSol,x_out,y_out);

    % built up the refined force field:
    refinedForceField(frame).pos=horzcat(x_out, y_out);
    refinedForceField(frame).vec=horzcat(fx, fy);

    refinedForceField(frame).par=forceField(frame).par;
    refinedForceField(frame).par.forceMesh = refinedForceMesh;
    refinedForceField(frame).par.sol_coef  = sol_coef;
    refinedForceField(frame).par.M         = M;


    figure(2)
    quiver(forceField(frame).pos(:,1),forceField(frame).pos(:,2),forceField(frame).vec(:,1)/maxForcePlot,forceField(frame).vec(:,2)/maxForcePlot,0,'r');
    hold on
    quiver(refinedForceField(frame).pos(:,1),refinedForceField(frame).pos(:,2),refinedForceField(frame).vec(:,1)/maxForcePlot,refinedForceField(frame).vec(:,2)/maxForcePlot,0,'b');
    hold off
    title(['Force field no: ',num2str(frame)])
    axis equal
    xlim([min_x max_x])
    ylim([min_y max_y])
    set(gca,'YDir','reverse')
    
end