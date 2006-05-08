function LG_plot(movieFrame, idlist, axesH, frameSizeMu, pixelSize, plotOptions)
%LG_plot plots data for labelgui2

% defaults
signList = 'o*h';
if ~isempty(idlist)
    colorMap = plotOptions.colorMap;
    popup = plotOptions.tagPopupMenuH;
end

% turn off warning of obsolete imshow syntax
warnstate = warning;
warning off Images:imshow:obsoleteSyntaxXY

% For less typing: We loop through the axesH (xy, yz, xz), generate the
% image and the positions. Then plot. This works well with additional
% windows.

for ix = 1:length(axesH)

    t = floor((ix-1)/3)+1;
    if ~isempty(idlist)
    nTags = size(idlist(t).linklist,1);
    end
    switch rem(ix-1,3) + 1
        case 1 % xy window

            % prepare data for plotting

            % image
            if isempty(plotOptions.currentSlice)
                % maximum projection
                plotImage = max(movieFrame(:,:,:,:,t),[],3);
            else
                % take current slice
                plotImage = movieFrame(:,:,plotOptions.currentSlice,:,t);
            end

            if ~isempty(idlist) && ~isempty(idlist(t).linklist)
                % points: x,y,spotNumber,flag. Keep positions in microns
                plotData = idlist(t).linklist(:,[9,10,2,3,5]);
            else
                plotData = [];
            end

            % remove old children
            currentAxes = axesH(ix);
            axesHandles = get(currentAxes);
            if ~isempty(axesHandles.Children)
                delete(axesHandles.Children);
            end
            
            axes(currentAxes);
            h = imshow([0,frameSizeMu(2)]+0.5*pixelSize(1),...
                [0,frameSizeMu(1)]+0.5*pixelSize(1),plotImage,[]);
            set(h,'Parent',currentAxes)
            axis image
            set(currentAxes,'NextPlot','add');
            

            

        case 2 % yz window

            % prepare data for plotting

            % plotImage
                % maximum projection
                plotImage = squeeze(max(movieFrame(:,:,:,:,t),[],2));

            if ~isempty(idlist) && ~isempty(idlist(t).linklist)
                % points: x,y,spotNumber,flag. Keep positions in microns
                plotData = idlist(t).linklist(:,[11,10,2,3,5]);
            else
                plotData = [];
            end

            % remove old children
            currentAxes = axesH(ix);
            axesHandles = get(currentAxes);
            if ~isempty(axesHandles.Children)
                delete(axesHandles.Children);
            end
            
            axes(currentAxes);
            h = imshow([0,frameSizeMu(3)]+0.5*pixelSize(2),...
                [0,frameSizeMu(1)]+0.5*pixelSize(1),plotImage,[]);
            set(h,'Parent',currentAxes)
            axis image
            set(currentAxes,'NextPlot','add');
            
          
        case 3 % xz window

            % prepare data for plotting

            % plotImage
                % maximum projection
                plotImage = squeeze(max(movieFrame(:,:,:,:,t),[],1))';

            if ~isempty(idlist) && ~isempty(idlist(t).linklist)
                % points: x,y,spotNumber,flag. Keep positions in microns
                plotData = idlist(t).linklist(:,[9,11,2,3,5]);
                
            else
                plotData = [];
            end

            % remove old children
            currentAxes = axesH(ix);
            axesHandles = get(currentAxes);
            if ~isempty(axesHandles.Children)
                delete(axesHandles.Children);
            end
            
            axes(currentAxes);
            h = imshow([0,frameSizeMu(2)]+0.5*pixelSize(1),...
                [0,frameSizeMu(3)]+0.5*pixelSize(2),plotImage,[]);
            set(h,'Parent',currentAxes)
            axis image
            %
    set(currentAxes,'NextPlot','add');

    end % switch

    % plot

    % show plotImage
    axes(currentAxes);
    
    % plot rectangles
    %
    %
    %
    
    % plotData is: ["x", "y", spotNumber, flag1, flag2]
    % prepare for fusions: mark the primary spot with flag -3
    
    

    % plot tags
    if ~isempty(plotData)
        
        fusionIdx = plotData(:,4) == 3;
    allFusionIdx = ismember(plotData(:,3),plotData(fusionIdx,3));
    plotData(allFusionIdx,4) = 3 * sign(fusionIdx(allFusionIdx)-0.5);
        
        for iTag = 1:nTags
            % if it's a genuine tag, plot first sign. If lost, plot second.
            % If single occurence, third
            switchNum = (plotData(iTag,3) > 0) + ...
                2*(plotData(iTag,5) == 2 || plotData(iTag,5) == 3) + ...
                4*(abs(plotData(iTag,4)) == 3);
                % 0: lost tag, not single occurence
                % 1: good spot, goodTag
                % 2: lost tag, single occurence
                % 3: good spot, single occurence
                % 5: fusion
            switch switchNum
                case 0 % lost tag, not from single occurence. Plot with parenthesis
                    ph = plot(currentAxes,...
                        plotData(iTag,1),plotData(iTag,2),...
                        signList(2),'MarkerSize',6);
                    th = text(...
                        'Position',...
                        [plotData(iTag,1)+0.2,plotData(iTag,2),0],...
                        'String',...
                        sprintf('%s',plotOptions.labelColor{iTag}),...
                        'Background',[0.2,0.2,0.2],...
                        'Parent',currentAxes);
                    set([th,ph],'uicontextmenu',popup,...
                        'Color',colorMap(iTag,:),'UserData',iTag);

                case 1 % good spot, good tag
                    ph = plot(currentAxes,...
                        plotData(iTag,1),plotData(iTag,2),...
                        signList(1),'MarkerSize',6);
                    th = text('Parent',currentAxes,...
                        'Background',[0.2,0.2,0.2],...
                        'Position',...
                        [plotData(iTag,1)+0.2,plotData(iTag,2),0],...
                        'String',...
                        sprintf('%s',plotOptions.labelColor{iTag}));
                    set([th,ph],'uicontextmenu',popup,...
                        'Color',colorMap(iTag,:),'UserData',iTag);

                case 2 % lost tag, single occurence. don't plot
                    ph = [];

                case 3 % single occurence
                    ph = plot(currentAxes,...
                        plotData(iTag,1),plotData(iTag,2),...
                        signList(3),'MarkerSize',6);
                    th = text(...
                        'Position',...
                        [plotData(iTag,1)+0.2,plotData(iTag,2),0],...
                        'Background',[0.2,0.2,0.2],...
                        'String',...
                        sprintf('%s',plotOptions.labelColor{iTag}),...
                        'Parent',currentAxes);
                    set([th,ph],'uicontextmenu',popup,...
                        'Color',colorMap(iTag,:),'UserData',iTag);
                case 5 % fusion. 
                    % sign of plotData(iTag,4) tells whether it's primary
                    % or secondary (the latter, if positive)
                    if sign(plotData(iTag,4)) > 0
                        % secondary fusion. No text
                        ph = plot(currentAxes,...
                        plotData(iTag,1),plotData(iTag,2),...
                        signList(1),'MarkerSize',8);
                        th = [];
                    else
                        ph = plot(currentAxes,...
                        plotData(iTag,1),plotData(iTag,2),...
                        signList(1),'MarkerSize',4);
                    th = text(...
                        'Position',...
                        [plotData(iTag,1)+0.2,plotData(iTag,2),0],...
                        'Background',[0.2,0.2,0.2],...
                        'String',...
                        'fusion',...
                        'Parent',currentAxes);
                    end
                    set([th,ph],'uicontextmenu',popup,...
                        'Color',colorMap(iTag,:),'UserData',iTag);
                    
                otherwise
                    h = errordlg('unknown plot option','error in LG_plot');
                    uiwait(h)
                    return
            end

        end
    end
end % loop

warning(warnstate)
