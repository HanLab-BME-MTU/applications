function [constrForceField]=calcElEnergies(constrForceField,forceField,frame,displField,meshPtsFwdSol)
% Solve for the displacement field of the forcefield and then
% calculate the elastic energy.
if nargin <5
    meshPtsFwdSol=[];
end

if ~isfield(constrForceField{frame},'cell')
    display(['No cell found for frame: ',num2str(frame),'!?'])
    return;
end
j=1;
while j<=length(constrForceField{frame}.cell)
    while isempty(constrForceField{frame}.cell{j}.pos) && j<=length(constrForceField{frame}.cell)
        display(['There is a bug that needs to be fixed, cell: ',num2str(j),' in frame: ',num2str(frame),' has no .pos-field!?']);
        display('Set the energy field to empty!?');
        constrForceField{frame}.cell{j}.stats.elEnergy= [];
        constrForceField{frame}.cell{j}.stats.ratioCCForceOverElEnergy= [];
        j=j+1;
    end
    % Here we need to calculte the displacements only at the positions where
    % the force field doesn't vanish. For the elastic energy we want to
    % calculate int[u(x)*f(x)] and the integrand is zero if f=0. Thus
    % determine the bounds of the force field (Round to pixel value):    
    xmin=min(constrForceField{frame}.cell{j}.boundary(:,1));
    xmax=max(constrForceField{frame}.cell{j}.boundary(:,1));
    ymin=min(constrForceField{frame}.cell{j}.boundary(:,2));
    ymax=max(constrForceField{frame}.cell{j}.boundary(:,2));
    
    % Interpolate the force field within the cell on every pixel:
    [pixPosMat_x,pixPosMat_y]=meshgrid(xmin:1:xmax,ymin:1:ymax);
    
    forcePix_x=griddata(forceField(frame).pos(:,1),forceField(frame).pos(:,2),forceField(frame).vec(:,1),pixPosMat_x,pixPosMat_y,'linear'); % should be linear
    forcePix_y=griddata(forceField(frame).pos(:,1),forceField(frame).pos(:,2),forceField(frame).vec(:,2),pixPosMat_x,pixPosMat_y,'linear'); % should be linear
    
    bwMask=constrForceField{frame}.cell{j}.mask(ymin:ymax,xmin:xmax);
        
    % .mask and .boundary should encode the very same pixel loaction, thus
    % the following should work:
    inForceMat_x=forcePix_x;
    inForceMat_y=forcePix_y;
    inForceMat_x(~logical(bwMask))=0;
    inForceMat_y(~logical(bwMask))=0;
    inForceVec=horzcat(inForceMat_x(:),inForceMat_y(:));

%     % To check that everything works fine:    
%     figure(123)
%     contour(pixPosMat_x,pixPosMat_y,abs(inUMat_y)>0,300)
%     hold on
%     plot(constrForceField{frame}.cell{j}.boundary(:,1),constrForceField{frame}.cell{j}.boundary(:,2),'sr');
%     plot(constrForceField{frame}.cell{j}.pos(:,1),constrForceField{frame}.cell{j}.pos(:,2),'dg');
%     hold off
    
    % Take this force field now as interpolant. Here we can use interp2
    % which is a bit faster than griddata (which in tunr peforms a bit
    % better than TriScatteredInterp):
    force_x=@(x,y) interp2(pixPosMat_x,pixPosMat_y,inForceMat_x,x,y,'*linear'); % should be linear
    force_y=@(x,y) interp2(pixPosMat_x,pixPosMat_y,inForceMat_y,x,y,'*linear'); % should be linear

    % calculate the displacement at exactly the same position where the
    % force field has been measured!
    % The result from fft and conv agree within a range less then 1%

    
    [ux_fine uy_fine x_fine y_fine meshPtsFwdSol]=fwdSolution([xmin xmax],[ymin ymax],constrForceField{frame}.par.yModu_Pa,xmin,xmax,ymin,ymax,force_x,force_y,'fft','noIntp',meshPtsFwdSol); %'conv'
    % The corresponding position are: constrForceField{frame}.cell{j}.pos;
    x0=constrForceField{frame}.cell{j}.pos(:,1);
    y0=constrForceField{frame}.cell{j}.pos(:,2);
    constrForceField{frame}.cell{j}.uFwdSol=horzcat(interp2(x_fine,y_fine,ux_fine,x0,y0,'*cubic'),interp2(x_fine,y_fine,uy_fine,x0,y0,'*cubic'));
    
    if isfield(constrForceField{frame}.cell{j}.stats,'method') && strcmp(constrForceField{frame}.cell{j}.stats.method,'noIntp')
        factor=constrForceField{frame}.par.gridSpacing^2*constrForceField{frame}.par.pixSize_mu^3/10^6;
        constrForceField{frame}.cell{j}.stats.elEnergy= 1/2*sum(sum(constrForceField{frame}.cell{j}.uFwdSol.*constrForceField{frame}.cell{j}.vec))*factor;
    else
        if ~isfield(constrForceField{frame}.cell{j}.stats,'method') || ~strcmp(constrForceField{frame}.cell{j}.stats.method,'pixIntp')
            display('Method is unclear, used "Intp" for calculating el. Energy');
            constrForceField{frame}.cell{j}.stats.method='unclear';
        end
        
        % Calculate the displacement Field at all pixel postions:
        iu_x=@(x,y) interp2(x_fine,y_fine,ux_fine,x,y,'*cubic'); % should be cubic
        iu_y=@(x,y) interp2(x_fine,y_fine,uy_fine,x,y,'*cubic');
        
        uPix_x=iu_x(pixPosMat_x,pixPosMat_y);
        uPix_y=iu_y(pixPosMat_x,pixPosMat_y);
        
        % Take only the part of the displacement field that is within the
        % cell mask:
        inUMat_x=uPix_x;
        inUMat_y=uPix_y;
        inUMat_x(logical(~bwMask))=0;
        inUMat_y(logical(~bwMask))=0;
        inUVec=horzcat(inUMat_x(:),inUMat_y(:));
                
        gridSpacing=1;
        factor=gridSpacing^2*constrForceField{frame}.par.pixSize_mu^3/10^6;
        constrForceField{frame}.cell{j}.stats.elEnergy= 1/2*sum(sum(inUVec.*inForceVec))*factor;
    end
    constrForceField{frame}.cell{j}.stats.ratioCCForceOverElEnergy=constrForceField{frame}.cell{j}.stats.resForce.mag/constrForceField{frame}.cell{j}.stats.elEnergy;        
    constrForceField{frame}.cell{j}.stats.meshPtsFwdSol=meshPtsFwdSol;        
    clear('ux','uy')
    
    % This part calculates an APPROXIMATE value for the elastic energy of
    % each cell assuming that the displacements measured in the footprint
    % of each cell are only due to forces exerted by the same cell. In
    % general, however, forces exerted by neighboring cells also give a
    % contribution to ALL measured displacements. The two results obtained
    % from above and below agree within a range of ~30%.
    if nargin>3
        iDisplCoordx = TriScatteredInterp(displField(frame).pos(:,1),displField(frame).pos(:,2),displField(frame).vec(:,1));
        iDisplCoordy = TriScatteredInterp(displField(frame).pos(:,1),displField(frame).pos(:,2),displField(frame).vec(:,2));

        constrForceField{frame}.cell{j}.iuExp(:,1) = iDisplCoordx(constrForceField{frame}.cell{j}.pos(:,1),constrForceField{frame}.cell{j}.pos(:,2));
        constrForceField{frame}.cell{j}.iuExp(:,2) = iDisplCoordy(constrForceField{frame}.cell{j}.pos(:,1),constrForceField{frame}.cell{j}.pos(:,2));
        constrForceField{frame}.par.filter         = displField(frame).par.filter;
        
        factor=constrForceField{frame}.par.gridSpacing^2*constrForceField{frame}.par.pixSize_mu^3/10^6;
        constrForceField{frame}.cell{j}.stats.elEnergyAprox= 1/2*sum(sum(constrForceField{frame}.cell{j}.iuExp.*constrForceField{frame}.cell{j}.vec))*factor;
        clear('ux','uy')
    else
        display('Attention: approximate values for the elEnergy are not calculated!')
    end
    j=j+1;
end

% function [Zout]=myGriddata(Xin,Yin,Zin,Xout,Yout,method)
%     Zout=griddata(Xin, Yin, Zin, Xout, Yout, method);
%     nanMat=isnan(Zout);
%     Zout(nanMat)=0;
    
%   This is even slower than griddata:    
%     Zfunction=TriScatteredInterp(Xin, Yin, Zin, method);
%     Zout=Zfunction(Xout, Yout);
%     nanMat=isnan(Zout);
%     Zout(nanMat)=0;
% end
end