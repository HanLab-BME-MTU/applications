function objective = trackerObjectiveFunction(movieFrame,sourceInfo,targetInfo,initialParameters,constants,finalParms,movieNoise)
% this is for testChromdynExperiment(1) only. We move the tags along
% dimension 1 to get a 2D output objective function

nTags = length(targetInfo);

% deltas
deltas = [-5:0.5:5];
nDeltas = length(deltas);

sigmaResidual2 = zeros(nDeltas,nDeltas,2);

for iDelta = 1:nDeltas
    for jDelta = 1:nDeltas
        parameters = initialParameters;
        parameters(1) = parameters(1) + deltas(iDelta);
        parameters(4) = parameters(4) + deltas(jDelta);

        %don't calculate gradient
        [sourceInfo, targetInfo] = ...
            extractIntensities(sourceInfo, targetInfo,...
            movieFrame, constants, parameters, 0);
        for iTag = 1:nTags
            residuals = ...
                targetInfo(iTag).deltaInt(targetInfo(iTag).goodIdx);
            dof = length(residuals)-length(parameters);
            % since we did a subtraction of two noisy signals, the
            % noise should increase by sqrt(2)
            % potential improvement: calculate robust second moment
            sigmaResidual2(iDelta,jDelta,iTag) = (sum(residuals.^2)/(dof*2));

            % f-test
            %                     fRatio = sigmaResidual2(iTag,iter+1)/movieNoise(currentTarget)^2;
            %                     fProb  = fcdf(fRatio,dof,numel(movieFrame)/2);
            %                     isSuccess(iTag) = fProb < 0.95;

        end
    end
end

figure,
subplot(2,2,1), contourf(deltas,deltas,sigmaResidual2(:,:,1)',50);
hold on, 
plot(0,0,'.r');
plot(finalParms(1)-initialParameters(1),...
    finalParms(4)-initialParameters(4),'.g');
colorbar('peer',gca)


subplot(2,2,2), contourf(deltas,deltas,sigmaResidual2(:,:,2)',50);
hold on, 
plot(0,0,'.r');
plot(finalParms(1)-initialParameters(1),...
    finalParms(4)-initialParameters(4),'.g');
colorbar('peer',gca)


ssr2 = sum(sigmaResidual2,3);
subplot(2,2,3), contourf(deltas,deltas,ssr2',50);
hold on, 
plot(0,0,'.r');
plot(finalParms(1)-initialParameters(1),...
    finalParms(4)-initialParameters(4),'.g');
colorbar('peer',gca)



% plot also success rate
fRatio = sigmaResidual2/movieNoise^2;
fprob = fcdf(fRatio,dof,numel(movieFrame)/2);
fprob = max(fprob,[],3); % the maximum probability decides about rejection
subplot(2,2,4),
[c,h]=contourf(deltas,deltas,fprob',20);
%set(h,'ShowText','off','TextStep',get(h,'LevelStep')*4)
hold on, 
plot(0,0,'.r');
plot(finalParms(1)-initialParameters(1),...
    finalParms(4)-initialParameters(4),'.g');
%colorbar('peer',gca)


objective = sigmaResidual2;