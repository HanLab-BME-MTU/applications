function [ h ] = GCAVisualsFilopodiaOverlaysFromFilterSet(filoInfo,imgSize,varargin)
%GCAVisualsFilopodiaOverlaysFromFilterSet
% INPUT:
% filoInfo
% imgSize
% filoFilterSet:  set of
% justExt
%%
% Notes: might always want to apply the filoFilterSet BEFORE Plotting
% Parameters/Optional are the plotValues, external vs no external,or plot
% by type,  if more than one type it will try to match the values to a type
% and plot at the base of that filo.
%% INPUTPARSER
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS
ip.addRequired('filoInfo');
ip.addRequired('imgSize');
ip.addParameter('filoFilterSet',[]); % note keep the filter set here as it makes it easier to plot original IDs if user desires
ip.addParameter('plotValues',[]);
ip.addParameter('justExt',1); % flag for internal filo versus ext filo plotting
ip.addParameter('colorEmbed','b');
ip.addParameter('colorFiloBranch','r');

ip.addParameteR('colorBranch','r'); % only applicable in branch mode
ip.addParameter('forceFiloText',0); % might be able to get rid of.
ip.addParameter('overlayName','filopodiaOverlay');
ip.addParameter('branchMode',false); % branch mode will by default color code by type
ip.addParameter('colorBranch','r');

%%
%set the defaults for filoFilterSet
if isempty(ip.Results.filoFilterSet);
    filoFilterSet = true(length(filoInfo),1);
else
    filoFilterSet = ip.Results.filoFilterSet;
end

if strcmpi(ip.Results.plotValues,'IDs')
    plotValues = 1:length(filoInfo); % plot the ID number
else
    plotValues = ip.Results.plotValues;
end

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
%
toAdd{1} = 'Ext_'; %
toAdd{2} = 'Int_';

plotTextAtBranches = false;

% Apply the filter
filoInfoFilt = filoInfo(filoFilterSet);
c = [ip.Results.colorEmbed;ip.Results.colorFiloBranch];

if ip.Results.branchMode == true % will use ip.Results.colorFiloBranch as the stem
    % find the smallest type
    types = vertcat(filoInfo(:).type);
    if length(unique(types))>2;
        display(['You are plotting in branch mode, but' ...
            'but the filter has more than two filoTypes: Check FilterSet']);
    end
    types(types==0) = 1; % treat any no branch filo as stem
    stemType = min(unique(types));
    filoInfoFilt= filoInfo(types==stemType); % set the filter info equal to the stem
    % plot the
    filoInfoBranch = filoInfo(types~=stemType);
    
    % test if the values equal the number of branches or the nubmer of
    % stems
    
    if length(filoInfoBranch) == length(plotValues);
        plotTextAtBranches = true;
    end
end % if ip.Results.branchMode == true
%% 
% Loop over the embedded/external filo if user desires.
for iType = typeStart:typeEnd
    colorC = c(iType,:);
    
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
        
        % plot the other pixels associated with the filo
        %         filoMask = zeros(imgSize);
        for i = 1:numel(filoInfoFilt)
            pixIndices = filoInfoFilt(i).([toAdd{iType} 'pixIndices']);
            idxEnd = find(pixIndices == filoInfoFilt(i).([toAdd{iType} 'endpointCoordFitPix']));
            pixIndicesPlot = pixIndices(1:idxEnd);
            [yC,xC] = ind2sub( imgSize  ,pixIndicesPlot);
            plot(xC,yC,'color',colorC,'Linewidth',2);
        end
        scatter(xycoordsEndFit(:,1),xycoordsEndFit(:,2),20,colorC,'filled');
        scatter(xycoordsBase(:,1),xycoordsBase(:,2),20, colorC,'filled');
        clear  idxEnd pixIndicesPlot xycoordsEndFit
    end % isempty
end    % iType 
%% Plot Branches if applicable 
% Ok the quickest way to do this right now is to repeat the above for branches - can maybe re-code if have time 
  %  the other option would be just to loop over the entire set with one
  % with one color and to only turn on the embed flag if user wants it . 
  % could just make below a subfunction 
   if ip.Results.branchMode == true 
       
    for i = 1:numel(filoInfoBranch)
            value = filoInfoFilt(i).('Ext_endpointCoordFitXY');
            if ~isnan(value);
                xycoordsEndFit(i,:) = filoInfoFilt(i).('Ext_endpointCoordFitXY');
            else
                xycoordsEndFit(i,:) = [NaN,NaN];
            end
            
        end
        % [yEndFit,xEndFit] = ind2sub(imgSize,coordsEndFitPix);
        % % xycoordsEndFit = [xEndFit yEndFit];
        xCoordsBase= arrayfun(@(x) x.('Ext_coordsXY')(1,1),filoInfoFilt);
        yCoordsBase = arrayfun(@(x) x.('Ext_coordsXY')(1,2),filoInfoFilt);
        xycoordsBase = [xCoordsBase' yCoordsBase'];
        % scatter(xycoordsEndFit(:,1),xycoordsEndFit(:,2),'r','filled');
        % scatter(xycoordsBase(:,1),xycoordsBase(:,2),'b','filled');
        %neuriteMask = analInfo(iFrame).masks.neuriteEdge;
        
        % plot the other pixels associated with the filo
        %         filoMask = zeros(imgSize);
        for i = 1:numel(filoInfoFilt)
            pixIndices = filoInfoFilt(i).( 'Ext_pixIndices');
            idxEnd = find(pixIndices == filoInfoFilt(i).( 'Ext_endpointCoordFitPix'));
            pixIndicesPlot = pixIndices(1:idxEnd);
            [yC,xC] = ind2sub( imgSize  ,pixIndicesPlot);
            plot(xC,yC,'color',ip.Results.colorBranch,'Linewidth',2);
        end
        scatter(xycoordsEndFit(:,1),xycoordsEndFit(:,2),20,ip.Results.colorBranch,'filled');
        scatter(xycoordsBase(:,1),xycoordsBase(:,2),20, ip.Results.colorBranch,'filled');
        clear  idxEnd pixIndicesPlot xycoordsEndFit
    end % branchMode    
       
   
      
%% Plot Text if applicable  
    if ~isempty(plotValues)
        if (length(plotValues) > 1 || forceFiloText == 1 )
            if plotTextAtBranches == 1
            arrayfun(@(i,j) text(i.Ext_coordsXY(1,1),i.Ext_coordsXY(1,2),numstr(j,3),'color','k'),filoInfoBranch,plotValues); 
            else 
            
            arrayfun(@(i,j) text(i.Ext_coordsXY(1,1),i.Ext_coordsXY(1,2),num2str(j,3),'color','k'),filoInfoFilt,plotValues);
            end 
        else 
            text(5,5,num2str(plotValues)); % single value parameter not per filopodia
        end
    end % isempty(plotValues) 
 % ip.Results.branchMode

 end 