% test the PDE model:

% Generate the mesh for the PDE problem (and overload p,e,t). First
% generate the boundary. That is four curves: curveL,curveT,curveR,curveB.
% The BC on curveL will be Dirichlet and Neumann on the remaining boundary.

% in general, first a closed curve from boundary to boundary has to be
% found that contains the interface of interest:
curveInterface=constrForceField{30}.interface{3}.pos;
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
j=1;
curveL{j}  =curveInterface;
curveTRB{j}=vertcat(piece3,piece1);
curveTRB{j}=removeDoublePoints(curveTRB{j});
% Make sure that the polygon is in counterclockwise order:
if ispolycw(vertcat(curveL{j}(:,1),curveTRB{j}(:,1)), vertcat(curveL{j}(:,2),curveTRB{j}(:,2)))
    curveL{j}  =flipud(curveL{j});
    curveTRB{j}=flipud(curveTRB{j});
end

% Boundary 2:
j=2;
curveL{j}  =flipud(curveInterface);
curveTRB{j}=piece2;
% Make sure that the polygon is in counterclockwise order:
if ispolycw(vertcat(curveL{j}(:,1),curveTRB{j}(:,1)), vertcat(curveL{j}(:,2),curveTRB{j}(:,2)))
    curveL{j}  =flipud(curveL{j});
    curveTRB{j}=flipud(curveTRB{j});
end

for j=1:2
% The part of the interface that falls into the dilated region should count
% to the stress free Boundary. Thus intersect this interface with the
% segmented cell boundary:
maskStruct=bwboundaries(maskSegm',4);
curveSegm =maskStruct{1}; % this should be the same as segmRes.curve!?
max_x=max(max(curveSegm(:,1)),max(curveL{j}(:,1)));
max_y=max(max(curveSegm(:,2)),max(curveL{j}(:,2)));
extMaskSegm=zeros(max_y,max_x);
indMask=sub2ind(size(extMaskSegm),curveSegm(:,2), curveSegm(:,1));
extMaskSegm(indMask)=1;    
[rowsMask,colsMask]=size(extMaskSegm);
figure(j)
imagesc(extMaskSegm)
    
indcurveL{j} = sub2ind(size(extMaskSegm), curveL{j}(:,2), curveL{j}(:,1));    
indHits{j}   = find(extMaskSegm(indcurveL{j}));

% Check if there are less then two points
if length(indHits{j})<2
    display('Sorry, something went wrong, please re-draw the interface')
    return;
else %There might be more than two points. Then take the ones with the 
     % largest distance:
    maxD{j}=0;
    for i=1:length(indHits{j})-1
        cand1{j}=curveL{j}(indHits{j}(i),:);
        cand2{j}=curveL{j}(indHits{j}(i+1),:);
        dc1c2{j}=sum((cand2{j}-cand1{j}).^2);        
        if dc1c2{j}>maxD{j} 
            innerXP1{j}=cand1{j};
            idInnerXP1{j}=indHits{j}(i);
            innerXP2{j}=cand2{j};
            idInnerXP2{j}=indHits{j}(i+1);
            maxD{j}=dc1c2{j};
        end
    end    
end

innerCurveL{j}=curveL{j}(idInnerXP1{j}:idInnerXP2{j},:);
extCurveTRB{j}=vertcat(curveL{j}(idInnerXP2{j}:end-1,:),vertcat(curveTRB{j},curveL{j}(2:idInnerXP1{j},:)));

% Now divide the curvesTRB into three pieces each (that is needed for
% myRectMesh):
dI{j}=floor(size(curveTRB{j},1)/3);
curveT{j}=extCurveTRB{j}(1:dI{j},:);
curveR{j}=extCurveTRB{j}(dI{j}:2*dI{j},:);
curveB{j}=extCurveTRB{j}(2*dI{j}:end,:);

figure(10+j)
plot(curveL{j}(:,1),curveL{j}(:,2),'.b')
hold on
plot(curveT{j}(:,1),curveT{j}(:,2),'or')
plot(curveR{j}(:,1),curveR{j}(:,2),'.g')
plot(curveB{j}(:,1),curveB{j}(:,2),'ok')
hold off

% Reduce the number of points along all curves. Take domain 1 first:
dSkip=10;
innerCurveL_short{j}=innerCurveL{j}(1:dSkip:end,:);
innerCurveL_short{j}(end,:)=innerCurveL{j}(end,:);
innerCurveL{j}=innerCurveL_short{j};

curveT_short{j}=curveT{j}(1:dSkip:end,:);
curveT_short{j}(end,:)=curveT{j}(end,:);
curveT{j}=curveT_short{j};

curveR_short{j}=curveR{j}(1:dSkip:end,:);
curveR_short{j}(end,:)=curveR{j}(end,:);
curveR{j}=curveR_short{j};

curveB_short{j}=curveB{j}(1:dSkip:end,:);
curveB_short{j}(end,:)=curveB{j}(end,:);
curveB{j}=curveB_short{j};

extCurveTRB_short{j}=extCurveTRB{j}(1:dSkip:end,:);
extCurveTRB_short{j}(end,:)=extCurveTRB{j}(end,:);
extCurveTRB{j}=extCurveTRB_short{j};


% Now we are ready to use myRectMesh:
[meshCom{j},femCom{j},femMat{j},p{j},e{j},t{j}]=myRectMesh(innerCurveL{j},curveT{j},curveR{j},curveB{j},0,0,1,6);

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
globYoung.val=((0*currMask)+10^7);%.*exp(-double(distMat)/lengthScale);


% Now calculate the solution;
u=assempde(b,p{j},e{j},t{j},c,a,f);

% output is interpolated to the node points p:
[u1{j},u2{j},exx{j},eyy{j},exy{j},sxx{j},syy{j},sxy{j},u1x{j},u2x{j},u1y{j},u2y{j}]=postProcSol(u,p{j},t{j});


% calculate the stress/forces exerted on the interface given as a curve 
% composed of line segments:
dSkipInt=2;
curveInterface_sparce{j}=innerCurveL{j}(1:dSkipInt:end,:);
[center{j},fx_sum{j},fy_sum{j},ftot_sum{j},fx_int{j},fy_int{j}]=calcIntfacialStress(curveInterface_sparce{j},sxx{j},syy{j},sxy{j},p{j},'nearest');
ftot_sum{j}
[1.0830e+04 2.1967e+03]*14^2
[1.1859e+04 1.6846e+03]*14^2

% calculate the stress on the imposed Dirichlet boundary:
[center_bd{j},fx_sum_bd{j},fy_sum_bd{j},ftot_sum_bd{j},fx_int_bd{j},fy_int_bd{j}]=calcIntfacialStress(extCurveTRB{j},sxx{j},syy{j},sxy{j},p{j},'nearest');
ftot_sum_bd
sum(constrForceField{30}.roi.vec)*14.^2

% Remove the global variables:
clear globForce


figure(30+j)
pdesurf(p{j},t{j},u1{j})

figure(40+j)
pdesurf(p{j},t{j},u2{j})

figure(50+j)
pdesurf(p{j},t{j},sxx{j})

figure(60+j)
pdesurf(p{j},t{j},syy{j})

figure(70+j)
pdesurf(p{j},t{j},sxy{j})

end

currentImage = double(imread('registered_2010_03_24_TFM_wellsB123C1_2_w4488 Band Pass_s25_t30.tif'));
dPix_Plot=50;
min_x=min(constrForceField{30}.roi.pos(:,1))-dPix_Plot;
max_x=max(constrForceField{30}.roi.pos(:,1))+dPix_Plot;
min_y=min(constrForceField{30}.roi.pos(:,2))-dPix_Plot;
max_y=max(constrForceField{30}.roi.pos(:,2))+dPix_Plot;
[rows,cols]=size(currentImage);
maxForcePlot=0.3*max(sqrt(fx_sum{j}.^2+fy_sum{j}.^2))/forceField(30).par.gridSpacing;

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


figure(90+j)
currentImage = double(imread('registered_2010_03_24_TFM_wellsB123C1_2_w4488 Band Pass_s25_t30.tif'));
imagesc(currentImage)
colormap('gray')
hold on
quiver(constrForceField{30}.roi.pos(:,1),constrForceField{30}.roi.pos(:,2),-constrForceField{30}.roi.vec(:,1),-constrForceField{30}.roi.vec(:,2),'-w');
plot(constrForceField{30}.interface{3}.pos(:,1),constrForceField{30}.interface{3}.pos(:,2),'-k');
plot(constrForceField{30}.segmRes.curve(:,1),constrForceField{30}.segmRes.curve(:,2),'k')
plot(constrForceField{30}.segmRes.curveDilated(:,1),constrForceField{30}.segmRes.curveDilated(:,2),'-b');
for j=1:2
    plot(center{j}(:,1),center{j}(:,2),'-r');
    quiver(center{j}(:,1),center{j}(:,2),fx_sum{j}/maxForcePlot,fy_sum{j}/maxForcePlot,0)
    quiver(center_bd{j}(:,1),center_bd{j}(:,2),fx_sum_bd{j}/maxForcePlot,fy_sum_bd{j}/maxForcePlot,0)
end
hold off
title(['Stress profile along interface, frame: ',num2str(30)])
axis equal
xlim([max([1 min_x]) min([cols,max_x])])
ylim([max([1 min_y]) min([rows,max_y])])
set(gca,'YDir','reverse')

figure(100+j)
imagesc(currentImage)
colormap('gray')
hold on
% quiver(constrForceField{30}.roi.pos(:,1),constrForceField{30}.roi.pos(:,2),-constrForceField{30}.roi.vec(:,1),-constrForceField{30}.roi.vec(:,2),'-w');
plot(constrForceField{30}.interface{3}.pos(:,1),constrForceField{30}.interface{3}.pos(:,2),'-k');
plot(constrForceField{30}.segmRes.curve(:,1),constrForceField{30}.segmRes.curve(:,2),'k')
plot(constrForceField{30}.segmRes.curveDilated(:,1),constrForceField{30}.segmRes.curveDilated(:,2),'-b');
for j=1:2
    quiver(p{j}(1,1:100:end),p{j}(2,1:100:end),u1{j}(1:100:end)',u2{j}(1:100:end)','-m');
    plot(center{j}(:,1),center{j}(:,2),'-r');
    quiver(center{j}(:,1),center{j}(:,2),fx_sum{j}/maxForcePlot,fy_sum{j}/maxForcePlot,0)
    quiver(center_bd{j}(:,1),center_bd{j}(:,2),fx_sum_bd{j}/maxForcePlot,fy_sum_bd{j}/maxForcePlot,0)
end
hold off
title(['Stress profile along interface, frame: ',num2str(30)])
axis equal
xlim([max([1 min_x]) min([cols,max_x])])
ylim([max([1 min_y]) min([rows,max_y])])
set(gca,'YDir','reverse')

