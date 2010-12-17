% test the PDE model:

% Generate the mesh for the PDE problem (and overload p,e,t). First
% generate the boundary. That is four curves: curveL,curveT,curveR,curveB.
% The BC on curveL will be Dirichlet and Neumann on the remaining boundary.

% in general, first a closed curve from boundary to boundary has to be
% found that contains the interface of interest:
curveInterface=constrForceField{30}.interface{2}.pos;
curveCellCluster=constrForceField{30}.segmRes.curveDilated;
maskSegm=constrForceField{30}.segmRes.mask;


% figure(1)
% plot(curveInterface(:,1),curveInterface(:,2),'.b')
% hold on
% plot(curveCellCluster(:,1),curveCellCluster(:,2),'or')
% hold off

% The interface cuts the cluster in two pieces, each will give us a stress
% profile. Get the boundaries for both pieces:

checkXpt1=0;
checkXpt2=0;
xpt1=curveInterface(1,:);
xpt2=curveInterface(end,:);
for i=1:size(curveCellCluster,1)
    pt=curveCellCluster(i,:);
    if xpt1(1)==pt(1) && xpt1(2)==pt(2)
        indXpt1=i;
        checkXpt1=checkXpt1+1;
    elseif xpt2(1)==pt(1) && xpt2(2)==pt(2)
        indXpt2=i;
        checkXpt2=checkXpt2+1;
    end
end
if checkXpt1~=1 || checkXpt2~=1
    display('Something went wrong')
    break
end

% Make sure that xpt1 is really the first crossing point:
if indXpt1>indXpt2
    curveInterface=flipud(curveInterface);
    indTemp=indXpt2;
    indXpt2=indXpt1;
    indXpt1=indXpt2;
    xpt1=curveInterface(1,:);
    xpt2=curveInterface(end,:);
end

% The cell cluster boundary is divided into three pieces:
piece1=curveCellCluster(      1:indXpt1,:);
piece2=curveCellCluster(indXpt1:indXpt2,:);
piece3=curveCellCluster(indXpt2:end ,:);


% Boundary 1:
curveL1  =curveInterface;
curveTRB1=vertcat(piece3,piece1);
curveTRB1=removeDoublePoints(curveTRB1);
% Make sure that the polygon is in counterclockwise order:
if ispolycw(vertcat(curveL1(:,1),curveTRB1(:,1)), vertcat(curveL1(:,2),curveTRB1(:,2)))
    curveL1  =flipud(curveL1);
    curveTRB1=flipud(curveTRB1);
end

% Boundary 2:
curveL2  =flipud(curveInterface);
curveTRB2=piece2;
% Make sure that the polygon is in counterclockwise order:
if ispolycw(vertcat(curveL2(:,1),curveTRB2(:,1)), vertcat(curveL2(:,2),curveTRB2(:,2)))
    curveL2  =flipud(curveL2);
    curveTRB2=flipud(curveTRB2);
end

% The part of the interface that falls into the dilated region should count
% to the stress free Boundary. Thus intersect this interface with the
% segmented cell boundary:
maskStruct=bwboundaries(maskSegm',4);
curveSegm =maskStruct{1}; % this should be the same as segmRes.curve!?
max_x=max(max(curveSegm(:,1)),max(curveL1(:,1)));
max_y=max(max(curveSegm(:,2)),max(curveL1(:,2)));
extMaskSegm=zeros(max_y,max_x);
indMask=sub2ind(size(extMaskSegm),curveSegm(:,2), curveSegm(:,1));
extMaskSegm(indMask)=1;    
[rowsMask,colsMask]=size(extMaskSegm);
imagesc(extMaskSegm)
    
indcurveL1 = sub2ind(size(extMaskSegm), curveL1(:,2), curveL1(:,1));    
indHits    = find(extMaskSegm(indcurveL1));

% Check if there are less then two points
if length(indHits)<2
    display('Sorry, something went wrong, please re-draw the interface')
    return;
else %There might be more than two points. Then take the ones with the 
     % largest distance:
    maxD=0;
    for i=1:length(indHits)-1
        cand1=curveL1(indHits(i),:);
        cand2=curveL1(indHits(i+1),:);
        dc1c2=sum((cand2-cand1).^2);        
        if dc1c2>maxD 
            innerXP1=cand1;
            idInnerXP1=indHits(i);
            innerXP2=cand2;
            idInnerXP2=indHits(i+1);
            maxD=dc1c2;
        end
    end    
end

innerCurveL1=curveL1(idInnerXP1:idInnerXP2,:);
extCurveTRB1=vertcat(curveL1(idInnerXP2:end-1,:),vertcat(curveTRB1,curveL1(2:idInnerXP1,:)));

% Now divide the curvesTRB into three pieces each (that is needed for
% myRectMesh):
dI1=floor(size(curveTRB1,1)/3);
curveT1=extCurveTRB1(1:dI1,:);
curveR1=extCurveTRB1(dI1:2*dI1,:);
curveB1=extCurveTRB1(2*dI1:end,:);

dI2=floor(size(curveTRB2,1)/3);
curveT2=curveTRB2(1:dI2,:);
curveR2=curveTRB2(dI2:2*dI2,:);
curveB2=curveTRB2(2*dI2:end,:);
% 
% figure(4)
% plot(curveL1(:,1),curveL1(:,2),'.b')
% hold on
% plot(curveT1(:,1),curveT1(:,2),'or')
% plot(curveR1(:,1),curveR1(:,2),'.g')
% plot(curveB1(:,1),curveB1(:,2),'ok')
% hold off


% Reduce the number of points along all curves. Take domain 1 first:
dSkip=10;
innerCurveL1_short=innerCurveL1(1:dSkip:end,:);
innerCurveL1_short(end,:)=innerCurveL1(end,:);
innerCurveL1=innerCurveL1_short;

curveT1_short=curveT1(1:dSkip:end,:);
curveT1_short(end,:)=curveT1(end,:);
curveT1=curveT1_short;

curveR1_short=curveR1(1:dSkip:end,:);
curveR1_short(end,:)=curveR1(end,:);
curveR1=curveR1_short;

curveB1_short=curveB1(1:dSkip:end,:);
curveB1_short(end,:)=curveB1(end,:);
curveB1=curveB1_short;

extCurveTRB1_short=extCurveTRB1(1:dSkip:end,:);
extCurveTRB1_short(end,:)=extCurveTRB1(end,:);
extCurveTRB1=extCurveTRB1_short;


% Now we are ready to use myRectMesh:
[meshCom,femCom,femMat,p,e,t]=myRectMesh(innerCurveL1,curveT1,curveR1,curveB1,0,0,1,6);

% Define the coefficient of the PDE model:
[b,c,a,f]=definePDEmodelMixed;


% Now define the Youngs modulus and the forces acting on the body:
global globForce
global globYoung
globForce.pos= forceField(30).pos;
globForce.vec=-forceField(30).vec;

% Within the footprint of the cell there is a high Young's modulus, outside
% of this area a lower Young's modulus has to be taken into account:
currMask=constrForceField{30}.segmRes.mask;
[rows, cols]=size(currMask);
[globYoung.xmat,globYoung.ymat] = meshgrid(1:cols,1:rows);
% Simple step function:
%globYoung.val=(999*currMask)+1;%currMask*0+10^6;
% % Exponential or linear decay:
lengthScale=1000;
distMat = bwdist(currMask);
globYoung.val=((0*currMask)+10^6);%.*exp(-double(distMat)/lengthScale);







% Now calculate the solution;
u=assempde(b,p,e,t,c,a,f);

% output is interpolated to the node points p:
[u1,u2,exx,eyy,exy,sxx,syy,sxy,u1x,u2x,u1y,u2y]=postProcSol(u,p,t);

% calculate the stress/forces exerted on the interface given as a curve 
% composed of line segments:
dSkipInt=2;
curveInterface_sparce=innerCurveL1(1:dSkipInt:end,:);
[center,fx_sum,fy_sum,ftot_sum]=calcIntfacialStress(curveInterface_sparce,sxx,syy,sxy,p,'nearest');
ftot_sum
[1.0830e+04 2.1967e+03]*14^2
[1.1859e+04 1.6846e+03]*14^2

% calculate the stress on the imposed Dirichlet boundary:
[center_bd,fx_sum_bd,fy_sum_bd,ftot_sum_bd,fx_int_bd,fy_int_bd]=calcIntfacialStress(extCurveTRB1,sxx,syy,sxy,p,'nearest');
ftot_sum_bd
sum(constrForceField{30}.roi.vec)*14.^2

% Remove the global variables:
clear globForce


figure(3)
pdesurf(p,t,u1)

figure(4)
pdesurf(p,t,u2)

figure(5)
pdesurf(p,t,sxx)

figure(6)
pdesurf(p,t,syy)

figure(7)
pdesurf(p,t,sxy)

currentImage = double(imread('registered_2010_03_24_TFM_wellsB123C1_2_w4488 Band Pass_s25_t30.tif'));
dPix_Plot=50;
min_x=min(constrForceField{30}.roi.pos(:,1))-dPix_Plot;
max_x=max(constrForceField{30}.roi.pos(:,1))+dPix_Plot;
min_y=min(constrForceField{30}.roi.pos(:,2))-dPix_Plot;
max_y=max(constrForceField{30}.roi.pos(:,2))+dPix_Plot;
[rows,cols]=size(currentImage);
maxForcePlot=0.3*max(sqrt(fx_sum.^2+fy_sum.^2))/forceField(30).par.gridSpacing;

figure(8)
imagesc(globYoung.val);
hold on
quiver(constrForceField{30}.roi.pos(:,1),constrForceField{30}.roi.pos(:,2),-constrForceField{30}.roi.vec(:,1),-constrForceField{30}.roi.vec(:,2),'-w');
plot(constrForceField{30}.segmRes.curve(:,1),constrForceField{30}.segmRes.curve(:,2),'k')
plot(constrForceField{30}.segmRes.curveDilated(:,1),constrForceField{30}.segmRes.curveDilated(:,2),'-b');
hold off
title(['Imposed Youngs modulus and inverted traction forces, frame: ',num2str(30)])
axis equal
xlim([max([1 min_x]) min([cols,max_x])])
ylim([max([1 min_y]) min([rows,max_y])])
set(gca,'YDir','reverse')


figure(9)
currentImage = double(imread('registered_2010_03_24_TFM_wellsB123C1_2_w4488 Band Pass_s25_t30.tif'));
imagesc(currentImage)
colormap('gray')
hold on
quiver(constrForceField{30}.roi.pos(:,1),constrForceField{30}.roi.pos(:,2),-constrForceField{30}.roi.vec(:,1),-constrForceField{30}.roi.vec(:,2),'-w');
plot(constrForceField{30}.interface{2}.pos(:,1),constrForceField{30}.interface{2}.pos(:,2),'-k');
plot(constrForceField{30}.segmRes.curve(:,1),constrForceField{30}.segmRes.curve(:,2),'k')
plot(constrForceField{30}.segmRes.curveDilated(:,1),constrForceField{30}.segmRes.curveDilated(:,2),'-b');
plot(center(:,1),center(:,2),'-r');
quiver(center(:,1),center(:,2),fx_sum/maxForcePlot,fy_sum/maxForcePlot,0)
quiver(center_bd(:,1),center_bd(:,2),fx_sum_bd/maxForcePlot,fy_sum_bd/maxForcePlot,0)
hold off
title(['Stress profile along interface, frame: ',num2str(30)])
axis equal
xlim([max([1 min_x]) min([cols,max_x])])
ylim([max([1 min_y]) min([rows,max_y])])
set(gca,'YDir','reverse')

figure(10)
imagesc(currentImage)
colormap('gray')
hold on
% quiver(constrForceField{30}.roi.pos(:,1),constrForceField{30}.roi.pos(:,2),-constrForceField{30}.roi.vec(:,1),-constrForceField{30}.roi.vec(:,2),'-w');
quiver(p(1,1:100:end),p(2,1:100:end),u1(1:100:end)',u2(1:100:end)','-m');
plot(constrForceField{30}.interface{2}.pos(:,1),constrForceField{30}.interface{2}.pos(:,2),'-k');
plot(constrForceField{30}.segmRes.curve(:,1),constrForceField{30}.segmRes.curve(:,2),'k')
plot(constrForceField{30}.segmRes.curveDilated(:,1),constrForceField{30}.segmRes.curveDilated(:,2),'-b');
plot(center(:,1),center(:,2),'-r');
quiver(center(:,1),center(:,2),fx_sum/maxForcePlot,fy_sum/maxForcePlot,0)
quiver(center_bd(:,1),center_bd(:,2),fx_sum_bd/maxForcePlot,fy_sum_bd/maxForcePlot,0)
hold off
title(['Stress profile along interface, frame: ',num2str(30)])
axis equal
xlim([max([1 min_x]) min([cols,max_x])])
ylim([max([1 min_y]) min([rows,max_y])])
set(gca,'YDir','reverse')