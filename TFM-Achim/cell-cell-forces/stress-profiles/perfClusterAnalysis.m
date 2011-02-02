function constrForceField=perfClusterAnalysis(constrForceField,forceField,frame)
constrForceField{frame}.clusterAnalysis=[];
% The pixelsize in mu for this frame is: 
pixSize_mu=constrForceField{frame}.par.pixSize_mu;

% This ensures that there the pde mesh is at least twice as dense as the
% force mesh:

% The pde grid should be at least 4 times denser than the force mesh (=2
% refinements);
expt   = 0;
frac  = 2^(-expt);
Hmax  = max(1,constrForceField{frame}.par.gridSpacing*frac);
numRef= max(0,round(2-expt));  

% Define the coefficient of the PDE model:
[b,c,a,f]=definePDEmodelAllNeumann;

% an alternative model is provided by:
%[b,c,a,f]=definePDEmodelAllDirichlet;

% Now generate the mesh for the PDE problem
curve_long=constrForceField{frame}.segmRes.curveDilated;

dPix=10;
bndCurve=curve_long(1:dPix:end,:);
% Sometimes it happens that first and last point are the same. In that case
% remove the last point of the bndCurve:
if compPts(bndCurve(1,:),bndCurve(end,:))
    bndCurve(end,:)=[];
end
[p,e,t,dl]=circMesh(bndCurve,numRef,Hmax); % for a denser mesh: circMesh(bndCurve,2);
% display('More than 1 refinements might be needed?!');

% Now define the Youngs modulus and the forces acting on the body:
global globForce
global globYoung
globForce.pos= forceField(frame).pos; % Since we take the measured forceField as interpolant it doesn't matter what the gridSpacing in constrForceField is!
globForce.vec=-forceField(frame).vec;

% Within the footprint of the cell there is a high Young's modulus, outside
% of this area a lower Young's modulus has to be taken into account:
currMask=constrForceField{frame}.segmRes.mask;
[rows, cols]=size(currMask);
[globYoung.xmat,globYoung.ymat] = meshgrid(1:cols,1:rows);
% Simple step function:
%globYoung.val=(999*currMask)+1;%currMask*0+10^6;
% % Exponential or linear decay:
lengthScale=10;
distMat = bwdist(currMask);
globYoung.val=((0*currMask)+10^7).*exp(-double(distMat)/lengthScale);

% Check if the force field covers the whole cell cluster:
idxConvH=convhull(forceField(frame).pos);
convH=forceField(frame).pos(idxConvH,:);
[pixConvH]=createPixelPath(convH);
% If the conv. hull cuts through the dilated cluster mask we run into
% problems:
maskDilated=constrForceField{frame}.segmRes.maskDilated;
linIdx = sub2ind(size(maskDilated), pixConvH(:,2), pixConvH(:,1));
errsum=sum(maskDilated(linIdx));


% Now calculate the solution;
tic;
u=assempde(b,p,e,t,c,a,f);
toc;

% Do a few rounds of adaptive refinement. This improves the solution for
% the cell-cell forces by less than 1% and is super slow! Is not worth it!
[u,p,e,t] = adaptmesh(dl,b,c,a,f,'Mesh',p,e,t,'Ngen',5);
toc;

% [u,p,e,t] = adaptmesh(dl,'boundaryCondition',pdePar{1},pdePar{2},pdePar{3},'Ngen',meshQuality,'Mesh',p,e,t,'Nonlin',nonLin,'Init',Ui);

% output is interpolated to the node points p:
[u1,u2,exx,eyy,exy,sxx,syy,sxy,u1x,u2x,u1y,u2y]=postProcSol(u,p,t);

% calculate the stress on the imposed Dirichlet boundary:
[center_bd,fx_sum_bd,fy_sum_bd,ftot_sum_bd,sx_sum_bd,sy_sum_bd,nVec_mean]=calcIntfacialStress(bndCurve,sxx,syy,sxy,p,pixSize_mu,'nearest');
    
% calculate the force at each interface, fill it into the network.
for j=1:length(constrForceField{frame}.network.edge)

    % calculate the stress/forces exerted on the interface given as a curve 
    % composed of line segments:
    curve_interface=constrForceField{frame}.network.edge{j}.intf;
    dPixIntf=curve_interface(1:dPix:end,:);
    [center,fx_sum,fy_sum,fc,sx_sum,sy_sum,nVec_mean,char]=calcIntfacialStress(dPixIntf,sxx,syy,sxy,p,pixSize_mu,'linear');
    
    [fc1,fc2]=assFcWithCells(constrForceField,frame,j,fc,char);

    % Don't save all values in constrForceField. That would become too
    % confusing! We only store what is necessary to easily use the functions
    % postProcSol and calcIntfacialStress.
    constrForceField{frame}.network.edge{j}.dPixIntf=dPixIntf; % the full length interface is found in .network.edge{j}.intf
    constrForceField{frame}.network.edge{j}.cntrs = center;   % Aufpunkte der Kraefte.
    constrForceField{frame}.network.edge{j}.f_vec = horzcat(fx_sum,fy_sum);
    constrForceField{frame}.network.edge{j}.s_vec = horzcat(sx_sum,sy_sum);
    constrForceField{frame}.network.edge{j}.fc1   = fc1; % belongs to nodes(1)
    constrForceField{frame}.network.edge{j}.fc2   = fc2; % belongs to nodes(2)
    constrForceField{frame}.network.edge{j}.fc    = fc; % this value is obsolate!
    constrForceField{frame}.network.edge{j}.n_Vec = nVec_mean;
    constrForceField{frame}.network.edge{j}.char  = char;
    constrForceField{frame}.network.edge{j}.errs  = errsum;
    
    
    % Do we need to store these values in here at all?
    constrForceField{frame}.clusterAnalysis.intf{j}.dPixIntf=dPixIntf; % the full length interface is found in .network.edge{j}.intf
    constrForceField{frame}.clusterAnalysis.intf{j}.cntrs  = center;   % Aufpunkte der Kraefte.
    constrForceField{frame}.clusterAnalysis.intf{j}.f_vec  = horzcat(fx_sum,fy_sum);
    constrForceField{frame}.clusterAnalysis.intf{j}.s_vec  = horzcat(sx_sum,sy_sum);
    constrForceField{frame}.clusterAnalysis.intf{j}.fc     = fc; % same as: network.edge{j}.fc
    constrForceField{frame}.clusterAnalysis.intf{j}.n_Vec  = nVec_mean; % this is the mean of all normal vectors on the interface
    constrForceField{frame}.clusterAnalysis.intf{j}.edgeNum=j;
end

constrForceField{frame}.clusterAnalysis.sol.u=u;    %These take too much space!
constrForceField{frame}.clusterAnalysis.mesh.p=p;   %These take too much space!
constrForceField{frame}.clusterAnalysis.mesh.e=e;   %These take too much space!
constrForceField{frame}.clusterAnalysis.mesh.t=t;   %These take too much space!
constrForceField{frame}.clusterAnalysis.coef.a=a;   
constrForceField{frame}.clusterAnalysis.coef.b=b;   
constrForceField{frame}.clusterAnalysis.coef.c=c;
constrForceField{frame}.clusterAnalysis.coef.f=f;
constrForceField{frame}.clusterAnalysis.errs  =errsum; % this checks if the force field covered the wohle cluster.
constrForceField{frame}.clusterAnalysis.par.dPix=dPix; 
%constrForceField{frame}.clusterAnalysis.par.globForce.pos=globForce.pos;
%constrForceField{frame}.clusterAnalysis.par.globForce.vec=globForce.vec;
constrForceField{frame}.clusterAnalysis.par.globYoung=globYoung;
%constrForceField{frame}.clusterAnalysis.par.globYoung.xmat=globYoung.xmat;  %These take too much space!
%constrForceField{frame}.clusterAnalysis.par.globYoung.ymat=globYoung.ymat;  %These take too much space!
%constrForceField{frame}.clusterAnalysis.par.globYoung.val=globYoung.val;
constrForceField{frame}.clusterAnalysis.bnd.bndCurve=bndCurve; % used to generate the mesh!

constrForceField{frame}.clusterAnalysis.bnd.pos=center_bd;
constrForceField{frame}.clusterAnalysis.bnd.f_vec = horzcat(fx_sum_bd,fy_sum_bd);
constrForceField{frame}.clusterAnalysis.bnd.s_vec = horzcat(sx_sum_bd,sy_sum_bd);
constrForceField{frame}.clusterAnalysis.bnd.f_tot = ftot_sum_bd;

% Remove the global variables:
clear globForce
clear globYoung
