function trackerObjectiveFunctionBig(movieFrame,sourceInfo,targetInfo,initialParameters,constants,finalParms,movieNoise)
% this is for getting the objective function around each tag of a two-tag
% movie. We only move one tag at a time, because 6D arrays get very big
% very quickly, and I'd like to get quarter-pixel resolution

nTags = length(targetInfo);

% deltas
deltas = -5:0.25:5;
nDeltas = length(deltas);

sigmaResidual2 = zeros(nDeltas,nDeltas,nDeltas,nTags);


for iTag = 1:nTags
    for iDelta = 1:nDeltas
        for jDelta = 1:nDeltas
            for kDelta = 1:nDeltas
                parameters = initialParameters;
                parameters(1 + (iTag-1)*3) = parameters(1 + (iTag-1)*3) + deltas(iDelta);
                parameters(2 + (iTag-1)*3) = parameters(2 + (iTag-1)*3) + deltas(jDelta);
                parameters(3 + (iTag-1)*3) = parameters(3 + (iTag-1)*3) + deltas(kDelta);

                %don't calculate gradient
                [sourceInfo, targetInfo] = ...
                    extractIntensities(sourceInfo, targetInfo,...
                    movieFrame, constants, parameters, 0);

                residuals = ...
                    targetInfo(iTag).deltaInt(targetInfo(iTag).goodIdx);
                dof = length(residuals)-length(parameters);
                % since we did a subtraction of two noisy signals, the
                % noise should increase by sqrt(2)
                % potential improvement: calculate robust second moment
                sigmaResidual2(iDelta,jDelta,kDelta,iTag) = (sum(residuals.^2)/(dof*2));

            end
        end
    end
end
% figure,
% subplot(2,2,1), contourf(deltas,deltas,sigmaResidual2(:,:,1)',50);
% hold on,
% plot(0,0,'.r');
% plot(finalParms(1)-initialParameters(1),...
%     finalParms(4)-initialParameters(4),'.g');
% colorbar('peer',gca)
%
%
% subplot(2,2,2), contourf(deltas,deltas,sigmaResidual2(:,:,2)',50);
% hold on,
% plot(0,0,'.r');
% plot(finalParms(1)-initialParameters(1),...
%     finalParms(4)-initialParameters(4),'.g');
% colorbar('peer',gca)
%
%
% ssr2 = sum(sigmaResidual2,3);
% subplot(2,2,3), contourf(deltas,deltas,ssr2',50);
% hold on,
% plot(0,0,'.r');
% plot(finalParms(1)-initialParameters(1),...
%     finalParms(4)-initialParameters(4),'.g');
% colorbar('peer',gca)



% plot also success rate
fRatio = sigmaResidual2/movieNoise^2;
fprob = fcdf(fRatio,dof,numel(movieFrame)/2);
% fprob = max(fprob,[],3); % the maximum probability decides about rejection
% subplot(2,2,4),
% [c,h]=contourf(deltas,deltas,fprob',20);
% set(h,'ShowText','on','TextStep',get(h,'LevelStep')*4)
% hold on,
% plot(0,0,'.r');
% plot(finalParms(1)-initialParameters(1),...
%     finalParms(4)-initialParameters(4),'.g');
% colorbar('peer',gca)
%
%
% objective = sigmaResidual2;

% assign in base

if ~evalin('base','exist(''ct'',''var'')')
    error('please assign the counter in the base workspace')
end
ct = evalin('base','ct');

assignin('base',sprintf('sigmaResidual2_%i',ct),sigmaResidual2);
assignin('base',sprintf('fProb_%i',ct),fprob);
assignin('base',sprintf('fRatio_%i',ct),fRatio);
assignin('base',sprintf('final_%i',ct),finalParms-initialParameters);
assignin('base',sprintf('noise_%i',ct),movieNoise);
assignin('base',sprintf('deltas',ct),deltas);
