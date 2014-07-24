function [traj]=randomWalk(numSteps,initPos,Dx,Dy,xDrift,yDrift)

randx=sign(rand(numSteps,1)-0.5);
randy=sign(rand(numSteps,1)-0.5);

% plot(0.5*tanh(0.1*((1:100)-25))+1)
% xDrift=xDrift*(0.5*(tanh(0.1*((1:numSteps)-35))+1))';

% if rand is larger 0.5 step right or up else left or down:
steps=[randx*Dx+xDrift,randy*Dy+yDrift];

traj=cumsum(vertcat(initPos,steps));