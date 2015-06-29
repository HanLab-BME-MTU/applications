function [ cost ] = costEndpointIntersection( E,I )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 
    row_i = 1;
    col_i = 2;
    label_i = 3;
    order_i = 4;
    orientation_i = 5;
    response_i = 6;
    
euclidean = hypot(E(row_i)-I(:,row_i),E(col_i)-I(:,col_i));
angle = atan((E(row_i)-I(:,row_i))./(E(col_i)-I(:,col_i)));

% angle should be restricted between -pi/2 and pi/2
% the angle difference should be the smaller of deltaAngle and pi -
% deltaAngle
deltaAngle = abs(I(:,orientation_i) - angle);
deltaAngle2 = abs(E(:,orientation_i) - angle);
deltaAngle = min([deltaAngle, pi - deltaAngle, deltaAngle2, pi - deltaAngle2]')';
% deltaAngle



cost = euclidean./cos(deltaAngle);

% cannot be part of the same connected component connected closer by 10
% pixels
cost(E(label_i) == I(:,label_i) & abs(E(order_i) - I(:,order_i)) < 10) = Inf;

% distance cutoff
% euclidean_T = 5*sqrt(2); 
euclidean_T = 7.0711;
cost(euclidean == 0) = Inf;
cost(euclidean > euclidean_T) = Inf;
% cost = euclidean;


end

