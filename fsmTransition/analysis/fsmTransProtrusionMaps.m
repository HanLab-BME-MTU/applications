function [transProtrusionMap2D,edgeProtrusionMap2D,indices]=fsmTransProtrusionMaps(transToEdgeDist,distancesTrans,distancesEdge,interpStep,indices)
% fsmTransProtrusionMaps creates 2D colorful maps of the protrusions calculated by fsmTransProtrusions
%
% SYNOPSIS    [transProtrusionMap2D,edgeProtrusionMap2D,indices]=
%                  fsmTransProtrusionMaps(transToEdgeDist,distancesTrans,distancesEdge,interpStep,indices)
%
% INPUT       transToEdgeDist : structure (t = 1:number of frames).y  : y coordinate of the knot of the transition spline at time t
%                                                                 .x  : x coordinate of the knot of the transition spline at time t
%                                                                 .dy : y distance at .y of the knot of the transition spline to 
%                                                                       a point on the edge spline at time t
%                                                                 .dx : x distance at .x of the knot of the transition spline to 
%                                                                       a point on the edge spline at time t
%             distancesTrans  : structure (t = 1:number of frames).y  : y coordinate of the knot of the spline at time t
%                                                                 .x  : x coordinate of the knot of the spline at time t
%                                                                 .dy : y displacement at .y of the knot of the spline from time t to t+1
%                                                                 .dx : x displacement at .x of the knot of the spline from time t to t+1                                     
%                               This is returned by fsmTransProtrusions.
%             distancesEdge   : same as distancesTrans for the edge.
%                               This is returned by fsmTransProtrusions.
%             interpStep      : [y x]: defines how many how dense the maps have to be interpolated in the spatial and temporal (see OUTPUT) 
%                                      dimensions; accepted values are {y,x} = 0 .. 1. interpStep = [0.5 0.5] means that the maps will be 
%                                      interpolated onto a 2x denser map (an interpolation point is introduced between two map points); 
%                                      interpStep = [0.2 0.2] results in a 5x denser map.
%
% OUTPUT      transProtrusionMap2D: a 2D map, with the number of frames (multiplied by 1/interStep(1)) on the vertical axis and the number of
%                                   knots (multiplied by 1/interStep(2)) on the horizontal axis.
%             edgeProtrusionMap2D : same as transProtrusionMap2D, for the edge.
%
% DEPENDENCES   fsmTransProtrusionMaps uses { }
%               fsmTransProtrusionMaps is used by { fsmTransition }
%
% REMARK     THIS FUNCTION ASSUMES THAT THE DIFFERENCE IN FRAMES DUE TO THE TIME-AVERAGING OF THE TRANSITION HAS ALREADY BEEN TAKEN INTO
%            ACCOUNT WHEN THE STRUCTURE OF SPLINES HAS BEEN GENERATED (this is the case if the function is called through the user 
%            interface FSMTRANSITION).
%
% Aaron Ponti, September 6th, 2004

if nargin==4
    indices=[1:length(distancesTrans(1).y)];
end

% Set time difference between transition and edge to 0
dTime=0;

% Debug flag
DEBUG=0

len=length(distancesTrans);

% Initialize maps
indices=unique(fix(indices)); % indices must be int
if max(indices)>length(distancesTrans(1).dy)
    indices=indices(1):length(distancesTrans(1).dy);
    fprintf(1,'Changed indices: indices=[%d:1:%d]\n',indices(1),indices(end));
end
transProtrusionMap2D=zeros(length(distancesTrans)-dTime,length(distancesTrans(1).dy(indices)));
edgeProtrusionMap2D=transProtrusionMap2D;

for i=1:len-dTime

    dyTall=distancesTrans(i).dy';          dyT=dyTall(indices);
    dxTall=distancesTrans(i).dx';          dxT=dxTall(indices);
    dyEall=distancesEdge(i+dTime).dy';     dyE=dyEall(indices);
    dxEall=distancesEdge(i+dTime).dx';     dxE=dxEall(indices);
    dyD=transToEdgeDist(i).dy';   
    dxD=transToEdgeDist(i).dx';   
    
    yTall=distancesTrans(i).y';         yT=yTall(indices);
    xTall=distancesTrans(i).x';         xT=xTall(indices);
    yEall=distancesEdge(i+dTime).y';    yE=yEall(indices);
    xEall=distancesEdge(i+dTime).x';    xE=xEall(indices);
    yD=transToEdgeDist(i).y'; 
    xD=transToEdgeDist(i).x';
    
    if DEBUG==1
        
        figure;
        plot(xTall,yTall,'m-');
        hold on
        plot(xT,yT,'r.');
        plot(xEall,yEall,'b-');
        plot(xE,yE,'k.'); 
        quiver(xT,yT,dxT,dyT,0);
        quiver(xE,yE,dxE,dyE,0);
        quiver(xD,yD,dxD,dyD,0);
        plot(xD+0.5*dxD,yD+0.5*dyD,'g.-','LineWidth',2);
        axis ij
        
    end
    
    % Find for every vector in T the closest in D
    DTD=createDistanceMatrix([yT' xT'],[yD' xD']);
    minTD=min(DTD,[],2);
    indxTD=zeros(length(yT),1);
    for j=1:length(yT)
        tmp=find(DTD(j,:)==minTD(j));
        indxTD(j)=tmp(1);
    end
    
    % Find for every vector in E the closest in D
    DED=createDistanceMatrix([yE' xE'],[yD' xD']);
    minED=min(DED,[],2);
    indxED=zeros(length(yE),1);
    for j=1:length(yE)
        tmp=find(DED(j,:)==minED(j));
        indxED(j)=tmp(1);
    end
    
    % Calculate distances for transition (t->t+1)
    dT=sqrt(dyT.^2+dxT.^2);
    
    % Calculate distances for edge (t->t+1)
    dE=sqrt(dyE.^2+dxE.^2);
    
    % Orient them with respect to 'transToEdge'
    sgnT=sign(dot([dyT;dxT],[dyD(indxTD);dxD(indxTD)],1));
    sgnT(find(sgnT==0))=1;
    dT=sgnT.*dT;
    %     sgnT=sign(dot([dyT;dxT],[-1*ones(1,length(dyT)); zeros(1,length(dyT))],1)); % Compare to vertical
    %     sgnT(find(sgnT==0))=1;
    %     dT=sgnT.*dT;
    
    sgnE=sign(dot([dyE;dxE],[dyD(indxED);dxD(indxED)],1));
    sgnE(find(sgnE==0))=1;
    dE=sgnE.*dE;
    %     sgnE=sign(dot([dyE;dxE],[-1*ones(1,length(dyE)); zeros(1,length(dyE))],1)); % Compare to vertical
    %     sgnE(find(sgnE==0))=1;
    %     dE=sgnE.*dE;
    
    % Store them into the 2D protrusion velocity map
    transProtrusionMap2D(i,1:length(dT))=dT;
    edgeProtrusionMap2D(i,1:length(dE))=dE;
    
end

% Interpolate
if interpStep~=0
    [x,y]=meshgrid(1:size(transProtrusionMap2D,2),1:size(transProtrusionMap2D,1));
    [xi,yi]=meshgrid(1:interpStep(2):size(transProtrusionMap2D,2),1:interpStep(1):size(transProtrusionMap2D,1));
    transProtrusionMap2D=interp2(x,y,transProtrusionMap2D,xi,yi);
    edgeProtrusionMap2D=interp2(x,y,edgeProtrusionMap2D,xi,yi);
end

mx=max([max(transProtrusionMap2D(:)) abs(min(transProtrusionMap2D(:))) max(edgeProtrusionMap2D(:)) abs(min(edgeProtrusionMap2D(:)))]);
figure;
imshow(transProtrusionMap2D,[-mx mx]); colormap('jet'); title('Transition');colorbar;xlabel('x');ylabel('time (frames)');
figure;
imshow(edgeProtrusionMap2D,[-mx mx]); colormap('jet'); title('Edge');colorbar;xlabel('x');ylabel('time (frames)');

if DEBUG==1
    
    figure;
    imshow(abs(edgeProtrusionMap2D)-abs(transProtrusionMap2D),[-mx mx]); colormap('jet'); title('|Edge|-|Trans|');colorbar;xlabel('x');ylabel('time (frames)');
    
    figure;surf(transProtrusionMap2D);colormap('jet');title('Trans');xlabel('x');ylabel('time (frames)');
    figure;surf(edgeProtrusionMap2D);colormap('jet');title('Edge');xlabel('x');ylabel('time (frames)');
    
    stdTransProtrusionMap2D=std(transProtrusionMap2D,[],1);
    stdEdgeProtrusionMap2D=std(edgeProtrusionMap2D,[],1);
    
    figure;
    plot(stdTransProtrusionMap2D,'k-');
    hold on
    plot(ones(1,size(transProtrusionMap2D,2))*mean(stdTransProtrusionMap2D),'b-');
    plot(stdEdgeProtrusionMap2D,'r-');
    plot(ones(1,size(edgeProtrusionMap2D,2))*mean(stdEdgeProtrusionMap2D),'m-');
    legend('std(x), trans','mean(std(x)), trans','std(x), edge','mean(std(x)), edge');
    
    fprintf(1,'Mean std(x), trans = %1.2f\n',mean(stdTransProtrusionMap2D));
    fprintf(1,'Mean std(x), edge  = %1.2f\n',mean(stdEdgeProtrusionMap2D));

end
