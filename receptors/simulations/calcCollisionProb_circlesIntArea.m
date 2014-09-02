function collisionProb = calcCollisionProb_circlesIntArea(pwDist,primaryNodeRadius,partnerNodeRadius)
%CALCCOLLISIONPROB_CIRCLESINTAREA is used to calculate a probability of
%collision for pairs that were determined not to have collided via the path
%based approach by receptorAggregationAlg_maxWeightedMatching_potColl.
%
%   INPUT:
%           pwDist:     pairwise distance between the center points of the
%                       paths taken by the pair of particles (nodes).
%
%           primaryNodeRadius:  the radius of one of the nodes taken as
%                               half the path length
%
%           partnerNodeRadius:  the radius of the second node taken as half
%                               the path length
%
%   OUTPUT:
%           collisionProb:  the calculated collision probability of nodes 1
%                           and 2. Circles for the two nodes with the radii
%                           given above can be formed and their
%                           intersection, if any, calculated using a
%                           formula (below). The collision probability is
%                           then calculated as the ratio of this area to
%                           the area of the larger of the two circles.
%
%   Robel Yirdaw, 08/27/14
%

    %Calculate respective areas covered by the primary and partner node.
    primaryNodeArea = pi*(primaryNodeRadius)^2;
    partnerNodeArea = pi*(partnerNodeRadius)^2;
    
    if ( ( (pwDist + primaryNodeRadius) < partnerNodeRadius ) ||...
            ( (pwDist + partnerNodeRadius) < primaryNodeRadius ) )
        %One circle is inside the other - taking the probability as the
        %ratio of the smaller circle to the larger.
        collisionProb = min(primaryNodeArea,partnerNodeArea) / max(primaryNodeArea,partnerNodeArea);
    else
    
        %The formula for the area formed by the intersectin of the two 
        %circles has three terms. 
        %(ref: mathworld.wolfram.com/Circle-CircleIntersection.html)
        
        areaTerm1 = (partnerNodeRadius^2)*acos((pwDist^2 + partnerNodeRadius^2 -...
            primaryNodeRadius^2)/(2*pwDist*partnerNodeRadius));

        areaTerm2 = (primaryNodeRadius^2)*acos((pwDist^2 + primaryNodeRadius^2 -...
            partnerNodeRadius^2)/(2*pwDist*primaryNodeRadius));

        areaTerm3 = 0.5*sqrt( (partnerNodeRadius + primaryNodeRadius - pwDist)*...
            (pwDist + partnerNodeRadius - primaryNodeRadius)*...
            (pwDist - partnerNodeRadius + primaryNodeRadius)*...
            (pwDist + partnerNodeRadius + primaryNodeRadius) );

        intersectionArea = areaTerm1 + areaTerm2 - areaTerm3;

        %The above calculation returns an imaginary value for the area
        %of intersection when pwDist > parimaryNodeRadius +
        %partnerNodeRadius. Extract real part.
        collisionProb = real(intersectionArea) / max(primaryNodeArea,partnerNodeArea);
        
    end
    
end
    
    