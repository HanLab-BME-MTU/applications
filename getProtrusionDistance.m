function protrusionDistance = getProtrusionDistance(protrusionVelocity, timeInterval)
sizeOfProtrusionVelocity = size(protrusionVelocity);
nWindows = sizeOfProtrusionVelocity(1);
nTimePoints = sizeOfProtrusionVelocity(2);
velocity=protrusionVelocity;

for j=1:nWindows
    distance=0;
    for k=1:nTimePoints
        if (isnan(velocity(j,k))==1)
            if (k==1)
                velocity(j,k)=0;
            else
                velocity(j,k)=velocity(j,k-1);
            end
        end
        distance = distance + velocity(j,k)*timeInterval;
        protrusionDistance(j,k)=distance;
    end
end
    