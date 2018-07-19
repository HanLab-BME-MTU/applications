function [neighbours]=getNeighbours(joint, skeletonImage)

% outputs x and y coordinates of the neighbour pixels in a joint

% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca


[rowUp, colUp]=find(skeletonImage(joint(1,2)-2:joint(1,2)-2, joint(1,1)-2:joint(1,1)+2)==1);
[rowLeft, colLeft]=find(skeletonImage(joint(1,2)-1:joint(1,2)+1, joint(1,1)-2:joint(1,1)-2)==1);
[rowDown, colDown]=find(skeletonImage(joint(1,2)+2:joint(1,2)+2, joint(1,1)-2:joint(1,1)+2)==1);
[rowRight, colRight]=find(skeletonImage(joint(1,2)-1:joint(1,2)+1, joint(1,1)+2:joint(1,1)+2)==1);

row=[rowUp'; rowLeft+1; rowDown'+4; rowRight+1];
col=[colUp'; colLeft; colDown'; colRight+4];


neighbours=[(col+joint(1,1)) (row+joint(1,2))]-3;

% Old version, when i was just removing the joint and not its vecinity
%[row, col]=find(skeletonImage(joint(1,2)-1:joint(1,2)+1, joint(1,1)-1:joint(1,1)+1)==1);
%neighbours=[col+joint(1,1) row+joint(1,2)]-2;

