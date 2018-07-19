function LG_plot(movieFrame, idlist, axesH, imageH, plotOptions)
%LG_plot plots data for labelgui2
%
% INPUT movieFrame : 3D array with intensities. Empty if the current and
%                    previous timepoint coincide (which means that we only
%                    need to replot the tags).
%       idlist     : current idlist
%       axesH      : handles to plot axes for tags
%       imageH     : handles to image axes
%           axes have been initialized so that they are in microns already
%       plotOptions: structure with fields
%                       maxTags = = movieWindowHandles.idlistData.maxTags;
%                       labelColor = movieWindowHandles.idlistData.labelcolor;
%                       colorMap = movieWindowHandles.colorMap;
%                       tagPopupMenuH = movieWindowHandles.tagPopupMenuH;
%                       currentSlice

%% defaults
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

%% read data

for ix = 1:length(axesH)

    t = floor((ix-1)/3)+1;
    if ~isempty(idlist)
        nTags = size(idlist(t).linklist,1);
    end
    switch rem(ix-1,3) + 1
        case 1 % xy window

            % prepare data for plotting

            % image
            if ~isempty(movieFrame)
                if isempty(plotOptions.currentSlice)
                    % maximum projection
                    plotImage = max(movieFrame(:,:,:,:,t),[],3);
                else
                    % take current slice
                    plotImage = movieFrame(:,:,plotOptions.currentSlice,:,t);
                end
            end

            if ~isempty(idlist) && ~isempty(idlist(t).linklist)
                % points: x,y,spotNumber,flag. Keep positions in microns
                plotData = idlist(t).linklist(:,[9,10,2,3,5]);
            else
                plotData = [];
            end


        case 2 % yz window

            % prepare data for plotting

            % plotImage
            % maximum projection
            if ~isempty(movieFrame)
                plotImage = permute(max(movieFrame(:,:,:,:,t),[],2),[1,3,2]);
            end

            if ~isempty(idlist) && ~isempty(idlist(t).linklist)
                % points: x,y,spotNumber,flag. Keep positions in microns
                plotData = idlist(t).linklist(:,[11,10,2,3,5]);
            else
                plotData = [];
            end


        case 3 % xz window

            % prepare data for plotting

            % plotImage
            % maximum projection
            if ~isempty(movieFrame)
                plotImage = permute(max(movieFrame(:,:,:,:,t),[],1),[2,3,1])';
            end

            if ~isempty(idlist) && ~isempty(idlist(t).linklist)
                % points: x,y,spotNumber,flag. Keep positions in microns
                plotData = idlist(t).linklist(:,[9,11,2,3,5]);

            else
                plotData = [];
            end

    end % switch

    %% plot

    % remove old children by removing overlay
    delete(findall(axesH(ix),'Tag','overlay'));

    % change CData of image to plot new image
    if ~isempty(movieFrame)
        set(imageH(ix),'CData',plotImage);
    end

    % show plotImage
    currentAxes = axesH(ix);
    %axes(currentAxes);

    % plot rectangles
    %
    %
    %

    % plotData is: ["x", "y", spotNumber, flag1, flag2]
    % prepare for fusions: mark the primary spot with flag -3



    % plot tags
    if ~isempty(plotData)



        for iTag = 1:nTags
            % if it's a genuine tag, plot first sign. If lost, plot second.
            % If single occurence, third
            switchNum = (plotData(iTag,3) > 0) + ...
                2*(plotData(iTag,5) == 2 || plotData(iTag,5) == 3) + ...
                4*(any(plotData(iTag,4) == [3,4]));
            % 0: lost tag, not single occurence
            % 1: good spot, goodTag
            % 2: lost tag, single occurence
            % 3: good spot, single occurence
            % 5: fusion
            switch switchNum
                case 0 % lost tag, not from single occurence. Plot with parenthesis
                    ph = plot(currentAxes,...
                        plotData(iTag,1),plotData(iTag,2),...
                        signList(2),'MarkerSize',6,...
                        'LineWidth',1);
                    th = text(...
                        'Position',...
                        [plotData(iTag,1)+0.2,plotData(iTag,2),0],...
                        'FontWeight','bold','String',...
                        sprintf('%s',plotOptions.labelColor{iTag}),...'Background',[0.2,0.2,0.2],...
                        'Parent',currentAxes);
                    set([th,ph],'uicontextmenu',popup,...
                        'Color',colorMap(iTag,:),'UserData',iTag,'Tag','overlay');

                case 1 % good spot, good tag
                    ph = plot(currentAxes,...
                        plotData(iTag,1),plotData(iTag,2),...
                        signList(1),'MarkerSize',6,...
                        'LineWidth',1);
                    th = text('Parent',currentAxes,...'Background',[0.2,0.2,0.2],...
                        'Position',...
                        [plotData(iTag,1)+0.2,plotData(iTag,2),0],...
                        'FontWeight','bold','String',...
                        sprintf('%s',plotOptions.labelColor{iTag}));
                    set([th,ph],'uicontextmenu',popup,...
                        'Color',colorMap(iTag,:),'UserData',iTag,'Tag','overlay');

                case 2 % lost tag, single occurence. don't plot
                    ph = [];

                case 3 % single occurence
                    ph = plot(currentAxes,...
                        plotData(iTag,1),plotData(iTag,2),...
                        signList(3),'MarkerSize',6,...
                        'LineWidth',1);
                    th = text(...
                        'Position',...
                        [plotData(iTag,1)+0.2,plotData(iTag,2),0],...%'Background',[0.2,0.2,0.2],...
                        'FontWeight','bold','String',...
                        sprintf('%s',plotOptions.labelColor{iTag}),...
                        'Parent',currentAxes);
                    set([th,ph],'uicontextmenu',popup,...
                        'Color',colorMap(iTag,:),'UserData',iTag,'Tag','overlay');
                case 5 % fusion.
                    % sign of plotData(iTag,4) tells whether it's primary
                    % or secondary (the latter, if positive)
                    if plotData(iTag,5) == 4
                        % secondary fusion. No text
                        ph = plot(currentAxes,...
                            plotData(iTag,1),plotData(iTag,2),...
                            signList(1),'MarkerSize',8,...
                            'LineWidth',1);
                        th = [];
                    else
                        ph = plot(currentAxes,...
                            plotData(iTag,1),plotData(iTag,2),...
                            signList(1),'MarkerSize',4,...
                            'LineWidth',1);
                        th = text(...
                            'Position',...
                            [plotData(iTag,1)+0.2,plotData(iTag,2),0],...'Background',[0.2,0.2,0.2],...
                            'FontWeight','bold','String',...
                            'fusion',...
                            'Parent',currentAxes,'Color',colorMap(iTag,:));
                    end
                    set([th,ph],'uicontextmenu',popup,...
                        'UserData',iTag,'Tag','overlay');

                otherwise
                    h = errordlg('unknown plot option','error in LG_plot');
                    uiwait(h)
                    return
            end

        end
    end
end % loop

warning(warnstate)
