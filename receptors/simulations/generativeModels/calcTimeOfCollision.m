function [timeOfCollision,distOfCollision] = calcTimeOfCollision(ptcle1,ptcle2)
%CALCTIMEOFCOLLISION uses initial and final positions of two particles to
%determine if particles have collided or approached each other. The time
%and distance of collision/closest approach is returned. The presence of a
%collision or closest approach is determined by looking at the translation
%of the particles in the x and y coordinates individually with respect to
%time.
%
%INPUT:     ptcle1/ptcle2: structures with initial and final x and y
%           coordinates of the two particles, given as fields xi,yi,xf,yf.
%
%OUTPUT:    timeOfCollision: if it exists, the time corresponding to the
%                            closest distance or collision of particles.
%                            Otherwise, NaN.
%
%           distOfCollision: if it exists, the minimum pairwise distance of
%                            the particles, which is 0 in case of collision
%                            or some nonzero value for closest approach.
%                            Otherwise, NaN.
%
%
%   Robel Yirdaw, 07/07/14
%

    %Initialize return values
    timeOfCollision = NaN;
    distOfCollision = NaN;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %Need to construct x(t) and y(t) for both particles. These are linear
    %functions. Get intercepts and slopes.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %The intercept values for both particles in x and y
    b1x = ptcle1.xi;
    b2x = ptcle2.xi;
    b1y = ptcle1.yi;
    b2y = ptcle2.yi;
    
    %Slopes for both particles in x and y - calculated in units of 1/dT,
    %where dt is the time step. Hence tf = 1 and ti = 0.
    m1x = (ptcle1.xf - ptcle1.xi)/1;
    m2x = (ptcle2.xf - ptcle2.xi)/1;
    m1y = (ptcle1.yf - ptcle1.yi)/1;
    m2y = (ptcle2.yf - ptcle2.yi)/1;
        
    %If needed:
    %Solve for individual time points in x and y where an intersection
    %occurs, i.e., the particles have the same coordinate in x and/or y.
    %tx = (b1x - b2x)/(m2x - m1x);
    %ty = (b1y - b2y)/(m2y - m1y);
    
    %During the time period, if the particles approach each other and then
    %move away, disance over time profile will have a minimum. If the
    %particles move away or towards each other for the entire time, then
    %the minimum pairwise distances will be at end or beginning of the
    %translation respectively. These are not considered here. Solve for
    %time of closest pairwise distance/collison:   
    tMin = ( (b2x - b1x)*(m1x - m2x) + (b2y - b1y)*(m1y - m2y) ) /...
        ( (m1x - m2x)^2 + (m1y - m2y)^2 );
    
    %If the minimum distance occurs during translation of the particles,
    %save it and also get the corresponding pairwise distance.
    if ( (tMin > 0) && (tMin < 1) )
        timeOfCollision = tMin;
        
        xDist = (m1x*tMin + b1x) - (m2x*tMin + b2x);
        yDist = (m1y*tMin + b1y) - (m2y*tMin + b2y);
        distOfCollision = sqrt(xDist^2 + yDist^2);
    end
        
    
    %{
    %To view the translations along x and y versus time:
    figure(); plot([0;1],[ptcle1.xi;ptcle1.xf],'k:sq',[0;1],[ptcle2.xi;ptcle2.xf],'r:o');
    figure(); plot([0;1],[ptcle1.yi;ptcle1.yf],'k--sq',[0;1],[ptcle2.yi;ptcle2.yf],'r--o');
    %}
end

