function probMerged = probMergedByDist(boxSize,diffConst,dT,distances,distancesPrior)


%% Main Body
distances = squeeze(distances);
diffConstStDev = sqrt(2*diffConst*dT);

%initialize probabilities
probMerged = zeros(length(distances),1);

%loop over all the provided distances and find the probability that, based
%only on distance, the particles are merged

noPriors = false;
if nargin < 5
    noPriors = true;
else
    distancesPrior = squeeze(distancesPrior);
end

for i = 1:length(distances)
    %the probability (density) that a separation distance of distance(i)
    %will be observed if the particles are merged, given the localization
    %error
    numMerged = normpdf(distances(i),0,diffConstStDev);
    
    if numMerged < 0.001
        probMerged(i) = 0;
        continue
    end
    
    %below is calculated the probability (density) that a separation
    %distance of distance(i) if the particles are not merged, given the
    %size of the box and possibly what the prior separation distance was
    if noPriors %then assume random placement within box
        if distances(i) <= boxSize
            numSplit = (-4*distances(i)/boxSize^3 + pi/boxSize^2 + distances(i)^2/boxSize^4).*2.*distances(i);
        else
            numSplit = -2/boxSize^2 + 4/boxSize^2*asin(boxSize/distances(i)) + 4/boxSize^3*sqrt(distances(i)^2-boxSize^2) - pi/boxSize^2 - distances(i).^2/boxSize^4;
            numSplit = numSplit*distances(i)*2;
        end
    else
        numSplit = 0;
        for theta = 0:pi/100:2*pi
            numSplit = numSplit + exp(-(distances(i)^2 + distancesPrior(i)^2 - 2*distances(i)*distancesPrior(i)*cos(theta))/2/diffConstStDev^2)*pi/100;
        end
        numSplit = numSplit*distances(i)/2/pi/diffConstStDev^2;
    end
    
    %this is the probability that two particles are merged based solely on
    %their separation.
    probMerged(i) = numMerged/(numSplit+numMerged);
end

%plot the points
%plot(probMerged(:,1),probMerged(:,2))
