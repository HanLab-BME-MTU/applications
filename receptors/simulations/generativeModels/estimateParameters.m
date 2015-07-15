function likelihood = estimateParameters(rTraj,gTraj,locError)


%% Parameters

dT = 0.01;
boxSize = 10;
var1 = locError^2;
var2 = 2*locError^2;
[numR,~,numFrames] = size(rTraj);
[numG,~,~] = size(gTraj);
kDiv = 40;

%% Body
likelihood = zeros(kDiv);

for r = 1:numR
    for g = 1:numG
        fprintf('\npair: %d-%d',r,g);
        separations = squeeze(sqrt(sum((rTraj(r,:,:)-gTraj(g,:,:)).^2,2)));
        
        P = zeros(2,2,numFrames);
        piProb = zeros(1,2);
        
        piProb(1) = separations(1)/var1*exp(-separations(1)^2/(2*var1));
        piProb(2) = 2*separations(1)/boxSize^2;
        
        %pCoeffLog = 0;
        for t = 2:numFrames
            if (separations(t)<8*sqrt(var1))
                P(1,1,t) = separations(t)/2/var1*exp(-separations(t)^2/(4*var1));
            end
            for theta = 0:pi/100:2*pi
                P(2,2,t) = P(2,2,t) + exp(-(separations(t)^2 + separations(t-1)^2 - 2*separations(t)*separations(t-1)*sin(theta))/2/var2)*pi/100;
            end
            P(2,2,t) = P(2,2,t)*separations(t)/2/pi/var2;
            if P(1,1,t) > 0
                if P(2,2,t) > 0
                    pMin = squeeze(sqrt(P(1,1,t)*P(2,2,t)));
                else
                    pMin = P(1,1,t);
                end
            else
                if P(2,2,t) > 0
                    pMin = P(2,2,t);
                else
                    pMin = realmin;
                    P(2,2,t) = realmin;
                end
            end
            P(:,:,t) = P(:,:,t)/pMin;
            %pCoeffLog = pCoeffLog + log(pMin);
        end
        
        
        
        for on = 1:kDiv
            for off = 1:kDiv
                pOn = on/kDiv;
                kOff = off/kDiv*2;
                
                T = zeros(2,2);  
                T(1,2) = kOff*dT;
                T(1,1) = 1 - T(1,2);
                
                L = [1 0; 0 1];
                for t = 2:numFrames
                    T(2,1) = (separations(t)<(8*sqrt(var1)))*pOn;
                    T(2,2) = 1 - T(2,1);
                    L = L*squeeze(P(:,:,t))*T;
                end
                
                L = piProb*L*ones(2,1);
                likelihood(on,off) = likelihood(on,off) + log(L);
            end
        end
    end
end

[~,maxIndx] = max(likelihood(:));
[bestPOn, bestKOff] = ind2sub(size(likelihood),maxIndx);

fprintf('\npOn = %g',bestPOn/kDiv);
fprintf('\nkOff = %g, pOff = %g',bestKOff/kDiv*2,bestKOff/kDiv*dT);