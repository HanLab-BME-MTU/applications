
function [coords] = getEndpoints(pixIdx,size,checkSingletons,getVector,endPointType)
% checkSingletons input added for maintaining functionality with
% skel2Graph4OutGrowthMetric - here need to just get coordinates of
%
% getVector: logical 1 or 0 if 1 this will get the local vector from
%           towards the endpoint: need to determine if
%
% endPointType: scalar 4 or 8 (CC)
%           An endpoint can be defined as pixels with only 1 neighboring CC
%           in a neighborhood of either 4 or 8 local pixels.

if nargin<3
    checkSingletons = 0;
end

if nargin<4
    getVector = 0;
end

if nargin<5
    endPointType = 8;
end

singleFlag = 0 ;

if checkSingletons == 1
    if length(pixIdx) == 1
        singleFlag = 1;
        
    end
    % singletons as well even though not technically an end point.
end

if singleFlag == 1
    
    % it is a singtlton
    [ye,xe]  = ind2sub(size,pixIdx);
    vect = [nan nan];
else
    
    
    maskC = zeros(size);
    maskC(pixIdx)=1;
    
    if endPointType == 8
        sumKernel = [1 1 1];
        % find endpoints of the floating candidates to attach (note in the
        % future might want to prune so that the closest end to the
        % body is the only one to be re-attatched: this will avoid double connections)
        endpoints = double((maskC.* (conv2(sumKernel, sumKernel', padarrayXT(maskC, [1 1]), 'valid')-1))==1);
        
    else
        endpoints =     bwmorph(maskC,'endpoints');
    end
    
    
    [ye,xe] = find(endpoints~=0);
    
    
    % need to remove
    if getVector == 1
        if ~isempty(ye) % sometimes there are weird structure with out end points in the NMS give these a vect of NaN
            for iPoint = 1:length(ye)
                dist = bwdistgeodesic(logical(maskC),xe(iPoint),ye(iPoint));
                % get the index of the 3rd pixel from the end point
                [y3Back,x3Back] = ind2sub(size,find(dist == 3));
                % calculate the vector in the direction toward the endpoint.
                if isempty(y3Back) % quick fix for now is to use a smaller piece back
                    [y3Back,x3Back] = ind2sub(size,find(dist==2)); % two pixels back
                    
                end
                
                vectX = (xe(iPoint)-x3Back);
                vectY = (ye(iPoint)-y3Back);
                %Added
                vectX = vectX(1) ;
                vectY = vectY(1);
                distC = sqrt(vectX^2 + vectY^2);
                vect(iPoint,1) = vectX/distC; % make it a unit vector
                vect(iPoint,2) = vectY/distC;
            end
        else
            vect = [NaN,NaN;NaN NaN];
            xe = [NaN ; NaN] ;
            ye = [NaN ;NaN] ;
        end % if ~isempty(ye)
    end
    
    sanityCheck = 0;
    if sanityCheck == 1
        figure
        imshow(maskC);
        
        hold on
        
        arrayfun(@(i) quiver(xe(i),ye(i),vect(i,1),vect(i,2),10),1:2);
    end
    
    
    
    
end
coords(:,1) = xe;
coords(:,2) = ye;
if getVector == 1;
    coords = [coords vect] ;% add the vector coordinates
end
end



