function trackerPlotBigObjectiveFunction
%TRACKERPLOTBIGOBJECTIVEFUNCTION plots the output of trackerObjectiveFunctionBig
%
% SYNOPSIS: trackerPlotBigObjectiveFunction
%
% INPUT The function will search for data in the workspace
%
% OUTPUT Only plots and direct writing to commandline so far
%
% REMARKS
%
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: Jonas Dorn
% DATE: 16-Feb-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get sigmaResidual2
sigmaResidual2Names=evalin('base','who(''regexp'',''sigmaResidual2*'')');

nObjectiveFunctions = length(sigmaResidual2Names);
if isempty(sigmaResidual2Names)
    error('no sigmaResidual2 found in base workspace')
end

% try to get other variables
deltas = evalin('base','deltas');
nDeltas = length(deltas);
fProbNames = evalin('base','who(''regexp'',''fProb*'')');
fRatioNames = evalin('base','who(''regexp'',''fRatio*'')');
finalNames = evalin('base','who(''regexp'',''final*'')');
noiseNames = evalin('base','who(''regexp'',''noise*'')');

% loop through sigmaResidual2 and plot.
% - nspx3 plot with x,y,z projections of the minimum of the objective
% function for each spot
% - histogram of the fRatio
% - disp min-coordinate, actual minimum, minfProb, fProb at minimum

for iO = 1:nObjectiveFunctions
    % get number of objective function
    sr2Name = sigmaResidual2Names{iO};
    objectiveIdx = regexp(sr2Name,'\d+$');
    objectiveNumber = sr2Name(objectiveIdx:end);
    
    % load sigmaResidual2
    sigmaResidual2 = evalin('base',sr2Name);
    % find finalDelta, fProb, fRatio
    finalDelta = evalin('base',finalNames{iO});
    fProb = evalin('base',fProbNames{iO});
    fRatio = evalin('base',fRatioNames{iO});
    noise = evalin('base',noiseNames{iO});
    
    
    if any(sigmaResidual2(:)<0)
        badIdx = sigmaResidual2<=0;
        sigmaResidual2(badIdx) = NaN;
        fProb(badIdx) = NaN;
        fRatio(badIdx) = NaN;
    end
    % find number of tags
    nTags = size(sigmaResidual2,4);
    
    % find extreme values of the objective function for the clim properties
    % of the axes
    maxRes = nanmax(sigmaResidual2(:));
    minRes = nanmin(sigmaResidual2(:));
    
    
    
    % launch figure and loop through tags, projections
    ofh = figure('Name',sprintf('optimization %s',objectiveNumber));
    for iTag = 1:nTags
        figure(ofh)
        currentFinal = finalDelta((iTag-1)*3+1:iTag*3);
        sri = sigmaResidual2(:,:,:,iTag);
        for iProj = 1:3
            % create axes
            ah = subplot(nTags,3,(iTag-1)*3+iProj);
            % set color limits
            set(ah,'CLim',[minRes,maxRes]);
            % make projection
            projection = squeeze(nanmin(sri,[],iProj));
            % plot
            contourf(ah,deltas,deltas,projection,20);
            % add start, end
            hold on
            finalXY = currentFinal;
            finalXY(iProj) = [];
            plot(0,0,'+g',finalXY(2),finalXY(1),'+r')
            % add colorbar
            colorbar('peer',ah)
        end % loop projection
        
        % calculate statistics
        
        % minimum coordinates
        [dummy,minIdx] = nanmin(sri(:));
        [mx,my,mz] = ind2sub([nDeltas,nDeltas,nDeltas],minIdx);
        minCoord = deltas([mx,my,mz]);
        fpi = fProb(:,:,:,iTag);
        [minProb, minProbIdx] = nanmin(fpi(:));
        [mx,my,mz] = ind2sub([nDeltas,nDeltas,nDeltas],minProbIdx);
        minProbCoord = deltas([mx,my,mz]);
        
        % display statistics
        disp(sprintf(['Obj %s, Tag %i\nMinCoord   %1.4f %1.4f %1.4f\n',...
            'FinalDelta %1.4f %1.4f %1.4f\n',...
            'MinProbPos %1.4f %1.4f %1.4f\n',...
            'MinProb %1.4f\nNoise %1.4f\n\n'], objectiveNumber, iTag, minCoord, ...
            currentFinal, minProbCoord, minProb,noise))
        
        % plot fRatios
%         figure('Name',sprintf('fRatios %s',objectiveNumber));
%         histogram(fRatio(:,:,:,iTag));
        
    end % loop tags
   
end
