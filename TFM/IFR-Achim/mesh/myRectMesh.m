function [meshCom,femCom,femMat,p,e,t]=myRectMesh(curveL,curveT,curveR,curveB,hIn,vIn,onlyMatlabFun,numRef)
% curves Nx2 Matrix: first column has to be x-coordinate, second column has
%                    to be the y-coordinate.
% points along curveL...curveB have to be in clockwise order!

% We have to define the curves as global variables to make them visible for
% the geometry M-file 'defineGeom'. This is ugly but the only work around.

%!!! Important: femMat is equivalent to e,p,t
%!!! But: meshCom and femCom are different from femMat!!!

%first check if input is correct:
[isClosed,curveL,curveT,curveR,curveB]=isClosedCurve(curveL,curveT,curveR,curveB);
if ~isClosed
    display('the given curves don t present a closed boundary')
    return;
end

global globalCurveL globalCurveT globalCurveR globalCurveB;

globalCurveL=curveL;
globalCurveT=curveT;
globalCurveR=curveR;
globalCurveB=curveB;

%**************************************************************************
% Use only Matlab functions to generate the mesh:
%**************************************************************************

% figure(1)
% pdegplot('geomDeformedRect'), axis equal

% figure(2)
[p,e,t]=initmesh('geomDeformedRect'); 
% pdemesh(p,e,t), axis equal

if onlyMatlabFun==1
    meshCom=[];
    femCom =[];
    femMat =[];
    if nargin==8 && ~isempty(numRef)
        for i=1:numRef
            [p,e,t]=refinemesh('geomDeformedRect',p,e,t);
        end
        figure(2)
        pdemesh(p,e,t), axis equal
    end
    return;
end

% Now create the mesh for Comsol using the comsol fuction meshenrich. This
% works:
femMat.mesh.p=p;
femMat.mesh.e=e;
femMat.mesh.t=t;
femMat.mesh = meshenrich(femMat.mesh);

%figure(3)
%meshplot(femMat);


%**************************************************************************
% Use Comsol functions to generate the mesh:
%**************************************************************************

femCom.geom=geomobject('geomDeformedRect');
%figure(5)
% This doesn't really work:
%geomplot(femCom,'Detail','fine'), axis equal;

% But the mesh is correct. Here we can easily implement irregular meshing
% along the different edges using the option 'hnumedg':
femCom.mesh = meshinit('geomDeformedRect','hnumedg',[1 2 3 4; vIn hIn vIn hIn]);

figure(6)
meshplot(femCom);
set(gca,'Ydir','reverse')

meshCom=femCom.mesh;

clear('globalCurveL','globalCurveT','globalCurveR','globalCurveB')


%**************************************************************************
% Here come the private functions:
%**************************************************************************

function [isClosed,curveL,curveT,curveR,curveB]=isClosedCurve(curveL,curveT,curveR,curveB)
    % Check if start and end points of the curves do agree!
    % we compare 2*4 coordinates, thus the sum has to be 8 if the curve is
    % closed:
    if sum((curveL(1:2,end)==curveT(1:2,1))+(curveT(1:2,end)==curveR(1:2,1))+(curveR(1:2,end)==curveB(1:2,1))+(curveB(1:2,end)==curveL(1:2,1)))==8
        % For this input the transposed has to be taken.
        isClosed=true;
        curveL=curveL';
        curveT=curveT';
        curveR=curveR';
        curveB=curveB';
    elseif sum((curveL(end,1:2)==curveT(1,1:2))+(curveT(end,1:2)==curveR(1,1:2))+(curveR(end,1:2)==curveB(1,1:2))+(curveB(end,1:2)==curveL(1,1:2)))==8
        isClosed=true;
    else
        isClosed=false;
    end
end

end



% To test the function use these lines of code:

% curveL=[0  -0
%         0  -1];
%   
% curveT=[0    -1
%         0.3  -1.1
%         0.5  -0.9
%         0.7  -0.6
%         0.9  -1.5
%         1    -1];
% 
% curveR=[1 -1
%         1 -0];
% 
% curveB=[1  -0
%         0  -0];
% 
%     
% [femCom,femMat,p,e,t]=myRectMesh(curveL,curveT,curveR,curveB,100,20);
