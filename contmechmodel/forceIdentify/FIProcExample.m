%This is an example script file that shows the steps for force identification
% and post display.

%Set up the project, the finit element modle and assemble the linear system to be solved.
setupForceProj;
calDispField;
setupModel;

%Solve the linear system.
solveLSBF;

%Post solution assemble and display.
postAssembleBF;
postDispAdhMap;
quiver(bfDisplayPx,bfDisplayPy,recDispU(:,1)*dispScale, ...
   recDispU(:,2)*dispScale,0,'y');
quiver(bfDisplayPx,bfDisplayPy,recBF(:,1)*bfScale, ...
   recBF(:,2)*bfScale,0,'r');
