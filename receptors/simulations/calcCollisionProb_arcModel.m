function collisionProb = calcCollisionProb_arcModel(rmsd,pwDist)
%CALCCOLLISIONPROB_ARCMODEL calculates the probability of two particles
%   colliding when their pairwise distance is greater than the association
%   distance but less than the corrected association distance. 
%
%   A particle can take a step size of rmsd on average. The possible steps
%   for a particle then make up a circle with radius = rmsd centered on the 
%   particle. The association distance is set by user while the corrected 
%   association distance is taken as 2*rmsd. Thus, when the pairwise distance 
%   is < 2*rmsd, the circles for the two particles will intersect. 
%   The probablity of collision is therefore determined by calculating the 
%   probaiblity of a particle landing along the arc inscribed inside the 
%   circle of the second particle and this is done for both 
%   particles - collisionProb = p1^2.
%
%       INPUT: 
%           rmsd:   root-mean-square-displacement for a particle
%           pwDist: pairwise distance of the two particles
%
%       OUTPUT:
%           collisionProb:  probability of collision
%
%   Robel Yirdaw, 06/04/14
%   
%
    %Determine angle between line connecting the two particles and radius
    %of circle to the intersection point
    theta = real(acosd((0.5*pwDist)/rmsd));
    %Get length of the arc formed by the angle above
    arcLength = (pi/180)*rmsd*theta;
    %Probability of collision as ratio of arc over circle's circumference,
    %for both particles
    collisionProb = ( (2*arcLength)/(2*pi*rmsd) )^2;
    
end