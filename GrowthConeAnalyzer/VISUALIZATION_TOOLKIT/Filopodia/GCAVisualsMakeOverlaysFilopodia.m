function [ h,filoInfoFilt ] = GCAVisualsMakeOverlaysFilopodia( filoInfo ,imgSize,filterByFit,justExt,c,plotOrientFlag,forceResponse)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% justExt is a flag for plotting only the external filopodia (looks less
% cluttered for presentation movies as I am trying to get this cleaned up)
switch justExt
    case 1
        typeEnd =1;
        typeStart =1;
    case 2
        typeEnd = 2;
        typeStart = 2;
    case 3
        typeStart  = 1;
        typeEnd = 2;
end


toAdd{1} = 'Ext_';
toAdd{2} = 'Int_';

if (nargin < 5 || isempty(c) )
    c  = ['r';'b'];
end 
if (nargin<6 || isempty(plotOrientFlag))
    
    plotOrientFlag = 0;
end

if (nargin<7 || isempty(forceResponse))
   forceResponse = 0;  
end


filoID = 1:length(filoInfo);
for iType = typeStart:typeEnd
    colorC = c(iType,:);
    if filterByFit == 1
        % filter out those that did not have good enough fits
        test = vertcat(filoInfo(:).([ toAdd{iType} 'exitFlag'])); % check exit flag - less than 1 means a majure failure
        % very infrequent though - can happen for a crossing check to see
        % what type of filopodia typical look like this.
        
        
        
        
        
        %logicals for filtering of filopodia
        toRemoveBasedOnExitFlag = (test<1 ); % bad fits
        
        
        % filter out any that might have passed the exitflag criteria but NOT
        % gave a number for the fit ==0  % maybe flag above later and any
        % internal filopodia that do not exist (these will likewise be
        % marked by an NaN in the datastructure)
        endpointCoordTest = vertcat(filoInfo(:).([ toAdd{iType} 'endpointCoordFitPix'])); % want this to filter out internal
        % put NaN for each filoInfo if not corresponding internal
        % counterpart.
        nanToRemove = isnan(endpointCoordTest);
        
        % Also remove based on the parameters of the fit ....typical those
        % with p(1) < 100 are crap (intensity without background)
        % however a user can define this. It is also sometimes helpful to
        % filter by the ratio p(1)/p(4) the intensity of the signal to the
        % intensity of the background- but for now I keep it first with just the
        % intensity above
        p1Values = arrayfun(@(x) x.([toAdd{iType} 'params'])(1),filoInfo);
        
        % sometimes I had a small problem where I did not make the .params
        % field equal to length 4 - here simply fix this so can extract
        % with the same ID - they should be filtered in nanToRemove anyway
        %toFix = arrayfun(@(x) isnan(x.([toAdd{iType} 'params'])(1)),filoInfo);
        %cellfun(@(x) toFix
        %p4Values = arrayfun(@(x) x.([toAdd{iType} 'params'])(4),filoInfo);
        %ratioP1toP4 = p1Values./p4Values;
        
        
        % FILTER BASED ON FITTED LENGTH
        lengthsPix = arrayfun(@(x) x.([toAdd{iType} 'length']),filoInfo);
       % toRemoveBasedOnP1 = p1Values<100;
       toRemoveBasedOnP1 = zeros(1,length(filoInfo)); 
        toRemoveBasedOnLength = lengthsPix<2;
        
        % quick fix for filo 3 
        %toRemoveBasedOnLength(3) =1;
        
        
        % FIND FILOPODIA THAT MIGHT NOT HAVE PASSED THE CRITERIA FOR AN
        % EFFECTIVE LENGTH METRIC BUT GAVE A LONG SIGNAL IN THE STEERABLE
        % FILTER 
        lengthsSteer = arrayfun(@(x) length(x.([toAdd{iType} 'pixIndicesBack'])),filoInfo);
        % plot in dashed lines if exit flag was < 1 or steer filter detect
        % > 9 pixels ~ 2 um
        idxBadButPlot = ((lengthsSteer>9)' & (toRemoveBasedOnP1' | toRemoveBasedOnExitFlag | toRemoveBasedOnLength'));
       % idxBadButPlot(3) = 1; % 
        
        % HAVE A DIFFERENT CRITERIA TO FILTER BRANCHES - THESE ARE MORE LIKELY 
        % NOISE AS YOU GO FURTHER OUT 
        type = arrayfun(@(x) x.type,filoInfo);
        idxBranches = type>1;
        %toRemoveCleanBranches = (idxBranches & lengthsSteer<4 & p1Values<300);
        toRemoveCleanBranches = zeros(1,length(filoInfo)); 
        
        % MAY ALSO WANT TO EVENTUALLY FILTER BY RESIDUALS.. 
        
        % FINAL DEFINITION OF INDICES TO KEEP
        idxKeep = ~ ( nanToRemove  | toRemoveBasedOnP1' | toRemoveBasedOnLength'| toRemoveCleanBranches'|toRemoveBasedOnExitFlag);
        
        % FILTER INTERNAL FILOPODIA WITH A LOW CONFIDENCE EXTERNAL COUNTERPART
        if strcmpi(toAdd{iType},'Int_'); % if interal filopodia consider
            % strength of the external candidate for determining
            % confidence
            % if the interal does not have a strong external counterpart do
            % not plot
            
            endpointCoordTest = vertcat(filoInfo(:).([ 'Ext_endpointCoordFitPix'])); % want this to filter out internal
            % put NaN for each filoInfo if not corresponding internal
            % counterpart.
            nanToRemove = isnan(endpointCoordTest);
            % also remove based on the parameters of the fit ....typical those
            % with p(1) < 1 are crap (intensity without background)
            % also filter by a rough estimate of signal to noise (ie p(1)/p(4))
            % ideally you want this to be > 1
            
            p1Values = arrayfun(@(x) x.(['Ext_params'])(1),filoInfo);
            %p4Values = arrayfun(@(x) x.([toAdd{iType} 'params'])(4),filoInfo);
            lengthsPix = arrayfun(@(x) x.(['Ext_length']),filoInfo);
            toRemoveBasedOnP1 = p1Values<100;
            toRemoveBasedOnLength = lengthsPix<2;
            idxKeepExt = ~ ( nanToRemove  | toRemoveBasedOnP1' | toRemoveBasedOnLength' );
            
            idxKeep = (idxKeepExt & idxKeep);
            
        end
        
        % GET THE STRUCTURE OF HIGH CONFIDENCE FILOPODIA
        % Plotted in bold 
        filoInfoFilt = filoInfo(idxKeep);
        
        % GET THE STRUCTURE OF LOW CONFIDENCE LENGTH MEASUREMENT BUT LIKELY REAL
        % SIGNAL - plotted as dashed 
        filoInfoOther = filoInfo(idxBadButPlot);
        
    else % plot everything in bold 
        filoInfoFilt  = filoInfo;
        if forceResponse == 0
        filoInfoOther = [];
        else 
            filoInfoOther = filoInfo; 
        
        end 
    end
    
    if ~isempty(filoInfoFilt)
        for i = 1:numel(filoInfoFilt)
            value = filoInfoFilt(i).([toAdd{iType} 'endpointCoordFitXY']);
            if ~isnan(value);
                xycoordsEndFit(i,:) = filoInfoFilt(i).([toAdd{iType} 'endpointCoordFitXY']);
            else
                xycoordsEndFit(i,:) = [NaN,NaN];
            end
            
        end
        % [yEndFit,xEndFit] = ind2sub(imgSize,coordsEndFitPix);
        % % xycoordsEndFit = [xEndFit yEndFit];
        xCoordsBase= arrayfun(@(x) x.([toAdd{iType} 'coordsXY'])(1,1),filoInfoFilt);
        yCoordsBase = arrayfun(@(x) x.([toAdd{iType} 'coordsXY'])(1,2),filoInfoFilt);
        xycoordsBase = [xCoordsBase' yCoordsBase'];
        % scatter(xycoordsEndFit(:,1),xycoordsEndFit(:,2),'r','filled');
        % scatter(xycoordsBase(:,1),xycoordsBase(:,2),'b','filled');
        %neuriteMask = analInfo(iFrame).masks.neuriteEdge;
        %roiYX = bwboundaries(neuriteMask);
        %cellfun(@(x) plot(x(:,2),x(:,1),'g'),roiYX);
        
        % plot the other pixels associated with the filo
        %         filoMask = zeros(imgSize);
        for i = 1:numel(filoInfoFilt)
            pixIndices = filoInfoFilt(i).([toAdd{iType} 'pixIndices']);
            idxEnd = find(pixIndices == filoInfoFilt(i).([toAdd{iType} 'endpointCoordFitPix']));
            pixIndicesPlot = pixIndices(1:idxEnd);
            [yC,xC] = ind2sub( imgSize  ,pixIndicesPlot);
            plot(xC,yC,'color',colorC,'Linewidth',1.5);
        end
        %             filoMask(pixIndicesPlot) = 1;
        %         end
        
        %    spy(filoMask,colorC)
        
        % plot the other pixels associated with the filo
        
        scatter(xycoordsEndFit(:,1),xycoordsEndFit(:,2),20,colorC,'filled');
        scatter(xycoordsBase(:,1),xycoordsBase(:,2),20, colorC,'filled');
        if plotOrientFlag == 1
            orient = vertcat(filoInfoFilt(:).orientation);
            % orient = arrayfun(@(x)
            %arrayfun(@(x) isempty(x),filoInfoFilt.orientation);
            text(xycoordsBase(:,1),xycoordsBase(:,2),num2str(orient,2),'color','y','FontSize',10);
            
        end
        % for i = 1:numel(filoInfoFilt)
        %     lengthMeas = filoInfoFilt(i).length;
        % text(xycoordsEndFit(i,1),xycoordsEndFit(i,2),num2str(lengthMeas,3),'color','y');
        % text(xycoordsBase(i,1),xycoordsBase(i,2),num2str(i),'color','g');
        % end
        
        % if isfield(filoInfoFilt,'conXYCoords');
        %   branches = vertcat(filoInfoFilt(:).conXYCoords);
        % if ~isempty(branches)
        %   scatter(branches(:,1),branches(:,2),20,'y','filled');
        % end
        % end
        
        %  % plot pixIndices used for fit
        %  fitPix = vertcat(filoInfo(:).pixIndicesUsedForFit) ;
        %  fitMask = zeros(imgSize);
        %  fitMask(fitPix) = 1;
        %  spy(fitMask,'y');
        
        
        clear  idxEnd pixIndicesPlot xycoordsEndFit
    end % isempty
    
    %plot the lower confidence filo
    if   ~isempty(filoInfoOther)
        [toPlotY,toPlotX]  =  arrayfun(@(x) ind2sub(imgSize,x.([toAdd{iType} 'pixIndicesBack'])),filoInfoOther,'uniformoutput',0);
        for i = 1:numel(toPlotY)
            plot(toPlotX{i}(:),toPlotY{i}(:),':','color',c(iType,:),'Linewidth',1.5);
        end
    end
end

