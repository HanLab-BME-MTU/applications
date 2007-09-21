function fuzzyGroupArma_visualize(data,typicality,doEllipse,centers, classification,dmCutoff,mdsScale,plt)
%FUZZYGROUPARMA_VISUALIZE is the utility for visualizing ARMA clusters
%
% Input: data (as in groupArmaDescriptors)
%        typicality (from fuzzy clustering)
%        doEllipse - do ellipse plot of AR,MA . (def:0)
%        centers (only useful if 2d)
%        classification: hard membership array (rows: data points. Cols:
%        clusters. 1 if data point i belongs to cluster j).
%        dmCutoff - value that should become the max distance (def: none)
%        mdsScale - scaling for the mds plot - def: none
%        plt - plot option. def: 1
%

% defaults

% test input
if nargin == 0 || isempty(data)
    data = groupArma_loadData;
end
if isempty(data)
    error('no data loaded')
else
    % data is data(1:n).fitResults.bla
    % move .bla up into data
    fn = fieldnames(data);
    nData = length(data);
    if any(strmatch(fn,'fitResults'))
        for iData = 1:nData
            fni = fieldnames(data(iData).fitResults);
            for f = 1:length(fni)
                data(iData).(fn{f}) = data(iData).fitResults(fn{f});
            end
            data(iData) = rmfield(data(iData),'fitResults');
        end
    end
end
if nargin < 2 || isempty(typicality)
    doTypicality = false;
else
    doTypicality = true;
end
if nargin < 3 || isempty(doEllipse)
    doEllipse = false;
end
if nargin < 4 || isempty(centers) || ~doEllipse
    centers = [];
end
if nargin < 5 || isempty(classification)
    classification = ones(nData,1);
end
if nargin < 6 || isempty(dmCutoff)
    dmCutoff = false;
end
if nargin < 7 || isempty(mdsScale)
    mdsScale = false;
end
if nargin < 8 || isempty(plt)
    plt = 1;
end

% make distance matrix for ctVAT
dm=zeros(nData);
for iData = 2:nData,
    for jData=1:iData-1
        dm(iData,jData) = armaxModelComp(data(iData),...
            data(jData));
        dm(jData,iData) = dm(iData,jData);
    end,
end

% make sure there are no zeros off-diagonal
dm(dm==0) = eps;
dm = dm .* (1-eye(nData));

% clusterTendencyVAT
[sortedIdx, dmSorted] = clusterTendencyVAT(dm);

% do nonlinear multidimensional scaling
try
    [scaledCoords, stress, disparities] = mdscale(dm, 2,'Criterion','stress');
    % calculate new dissimilarities for comparison
    scaledDm = pdist(scaledCoords);
    scaledDm = squareform(scaledDm);
    useLin = false;
    crit = 's';
catch
    try
        [scaledCoords, stress, disparities] = mdscale(dm, 2,'Criterion','sstress');
        % calculate new dissimilarities for comparison
        scaledDm = pdist(scaledCoords);
        scaledDm = squareform(scaledDm);
        crit = 'ss';
        useLin = false;
    catch
        useLin = true;
    end
end

% sammon mds
[scaledCoords2, stress2, disparities2] = mdscale(dm, 2,'Criterion','sammon');
% calculate new dissimilarities for comparison
scaledDm2 = pdist(scaledCoords2);
scaledDm2 = squareform(scaledDm2);

if useLin
    scaledCoords = scaledCoords2;
    stress = stress2;
    scaledDm = scaledDm2;
    disparities = disparities2;
    crit = 'lin';
end

% if doEllipse: get "ARMA-coords"
% do anyway, in order to properly align data
try
    % collect data
    arma = [catStruct(1,'data.arParamK(1)'),catStruct(1,'data.maParamK(1)')];
    covArma = cat(3,data.varCovMatF);

    % align mds. Do not scale, so that maximum distance is preserved
    [dummy,dummy,transform] = procrustes(arma,scaledCoords);
    scaledCoords = scaledCoords*transform.T;
    [dummy,dummy,transform2] = procrustes(arma,scaledCoords2);
    scaledCoords2 = scaledCoords2*transform2.T;

catch
    disp('fGA_visualize: no ellipse/alignment possible %s',lasterr)
    doEllipse = false;
end

% transform mds coords so that they start from 0.5
scaledCoords(:,1) = scaledCoords(:,1) - min(scaledCoords(:,1)) + 0.5;
scaledCoords(:,2) = scaledCoords(:,2) - min(scaledCoords(:,2)) + 0.5;




% find max distance
if dmCutoff
    dmScale = dmCutoff;
else
    dmScale = max(dm(:));
end
% find max mds
if mdsScale
    axMax = mdsScale + 0.5;
else
    axMax = max(scaledCoords(:))+0.5;
end


switch plt
    case 1
        % plot all
        figure('Name','VAT and MDS')
        ah1= subplot(2,2,1);
        imshow(dm,[0,dmScale]);
        set(ah1,'xtick',[],'ytick',[],'Visible','on','TickDir','in');
        set(get(ah1,'Title'),'String','unsorted DM');



        ah2 = subplot(2,2,3);
        imshow(dmSorted,[0,dmScale]);
        % create shifted y-labels
        yTickLabel = '';
        for i=nData:-1:1
            if isEven(i)
                yTickLabel(i,:) = sprintf('     %3.0f',sortedIdx(i));
            else
                yTickLabel(i,:) = sprintf('%3.0f     ',sortedIdx(i));
            end
        end

        set(ah2,'xtick',1:nData,'xtickLabel',[],...
            'ytick',1:nData,'ytickLabel',yTickLabel,'Visible','on','TickDir','in');
        %set(get(ah2,'Title'),'String','sorted DM');
        set(get(ah2,'Title'),'String',sprintf('sorted DM %2.3f - %2.3f',min(dm(:)),max(dm(:))));

        ah3 = subplot(2,2,2);% mds plot with points centered
        cMax = max(scaledCoords,[],1);
        xShift = ((axMax - 0.5) - cMax(1))/2;
        yShift = ((axMax - 0.5) - cMax(2))/2;
        plot(scaledCoords(:,1)+xShift,scaledCoords(:,2)+yShift,'.')
        text(scaledCoords(:,1)+xShift+0.15,scaledCoords(:,2)+ yShift,num2str((1:nData)'));
        set(get(ah3,'Title'),'String',sprintf('nonlinear MDS (%s)',crit));
        axis equal
        xlim([0,axMax]);
        ylim([0,axMax]);

        % shepard plot
        ah4 = subplot(2,2,4);
        [dummy,sortIdx] = sort(dm(:));
        plot(dm(:),scaledDm(:),'.',dm(sortIdx),disparities(sortIdx),'-r')
        set(get(ah4,'Title'),'String',sprintf('Stress: %1.3f',stress));
        xlabel('dissimilarities')
        ylabel('transformed distances')
        axis square

        if doEllipse
            % plot second figure
            figure('Name','Ellipse & sammon-mds')
            ah1 = subplot(2,3,[1,2,4,5]);
            plot(arma(:,1),arma(:,2),'+');
            xlabel('1st AR')
            ylabel('1st MA')
            text(arma(:,1),arma(:,2),num2str((1:nData)'));
            hold on
            for i=1:nData
                errorEllipse(covArma(:,:,i),arma(i,:),'conf',-1);
            end

            ah3 = subplot(2,3,3);
            cMax = max(scaledCoords2,[],1);
            xShift = ((axMax - 0.5) - cMax(1))/2;
            yShift = ((axMax - 0.5) - cMax(2))/2;
            plot(scaledCoords2(:,1)+xShift,scaledCoords2(:,2)+yShift,'.')
            text(scaledCoords2(:,1)+xShift+0.15,scaledCoords2(:,2)+ yShift,num2str((1:nData)'));
            set(get(ah3,'Title'),'String','MDS');
            axis equal
            xlim([0,axMax]);
            ylim([0,axMax]);




            ah4 = subplot(2,3,6);
            [dummy,sortIdx] = sort(dm(:));
            plot(dm(:),scaledDm2(:),'.',dm(sortIdx),disparities2(sortIdx),'-r')
            set(get(ah4,'Title'),'String',sprintf('Stress: %1.3f',stress2));
            xlabel('dissimilarities')
            ylabel('transformed distances')
            axis square
        end

        % for paper

    case 2
        % plot all
        figure


        ah1 = gca;
        plot(arma(:,1),arma(:,2),'+');
        xlabel('1st AR')
        ylabel('1st MA')
        text(arma(:,1),arma(:,2),num2str((1:nData)'));
        hold on
        for i=1:nData
            errorEllipse(covArma(:,:,i),arma(i,:),'conf',-1);
        end
        axis square

        figure
        ah2 = gca;
        imshow(dmSorted,[0,dmScale]);
        % create shifted y-labels
        yTickLabel = '';
        for i=nData:-1:1
            if isEven(i)
                yTickLabel(i,:) = sprintf('     %3.0f',sortedIdx(i));
            else
                yTickLabel(i,:) = sprintf('%3.0f     ',sortedIdx(i));
            end
        end

        set(ah2,'xtick',1:nData,'xtickLabel',[],...
            'ytick',1:nData,'ytickLabel',yTickLabel,'Visible','on','TickDir','in','Position',[0.1,0.1,0.8,0.8]);
        %set(get(ah2,'Title'),'String','sorted DM');
        set(get(ah2,'Title'),'String',sprintf('sorted DM %2.3f - %2.3f',min(dm(:)),max(dm(:))));

        figure % mds plot with points centered
        ah3 = gca;
        cMax = max(scaledCoords,[],1);
        xShift = ((axMax - 0.5) - cMax(1))/2;
        yShift = ((axMax - 0.5) - cMax(2))/2;
        plot(scaledCoords(:,1)+xShift,scaledCoords(:,2)+yShift,'.')
        text(scaledCoords(:,1)+xShift+0.15,scaledCoords(:,2)+ yShift,num2str((1:nData)'));
        set(get(ah3,'Title'),'String',sprintf('nonlinear MDS (%s)',crit));
        axis equal
        xlim([0,axMax]);
        ylim([0,axMax]);

        % shepard plot
        figure,
        ah4 = gca;
        [dummy,sortIdx] = sort(dm(:));
        plot(dm(:),scaledDm(:),'.',dm(sortIdx),disparities(sortIdx),'-r')
        set(get(ah4,'Title'),'String',sprintf('Stress: %1.3f',stress));
        xlabel('dissimilarities')
        ylabel('transformed distances')
        axis square


    case 3
        % plot 2
        figure
        ah2=gca;
        imshow(dmSorted,[0,dmScale]);
        % create shifted y-labels
        yTickLabel = '';
        for i=nData:-1:1
            if isEven(i)
                yTickLabel(i,:) = sprintf('     %3.0f',sortedIdx(i));
            else
                yTickLabel(i,:) = sprintf('%3.0f     ',sortedIdx(i));
            end
        end

        set(ah2,'xtick',1:nData,'xtickLabel',[],...
            'ytick',1:nData,'ytickLabel',yTickLabel,'Visible','on','TickDir','in','Position',[0.1,0.1,0.8,0.8]);
        %set(get(ah2,'Title'),'String','sorted DM');
        set(get(ah2,'Title'),'String',sprintf('sorted DM %2.3f - %2.3f',min(dm(:)),max(dm(:))));

        figure % mds plot with points centered
        ah3 = gca;
        cMax = max(scaledCoords,[],1);
        xShift = ((axMax - 0.5) - cMax(1))/2;
        yShift = ((axMax - 0.5) - cMax(2))/2;
        plot(scaledCoords(:,1)+xShift,scaledCoords(:,2)+yShift,'.')
        text(scaledCoords(:,1)+xShift+0.15,scaledCoords(:,2)+ yShift,num2str((1:nData)'));
        set(get(ah3,'Title'),'String',sprintf('nonlinear MDS (%s)',crit));
        axis equal
        xlim([0,axMax]);
        ylim([0,axMax]);

        % shepard plot
        figure
        ah4 = gca;
        [dummy,sortIdx] = sort(dm(:));
        plot(dm(:),scaledDm(:),'.',dm(sortIdx),disparities(sortIdx),'-r')
        set(get(ah4,'Title'),'String',sprintf('Stress: %1.3f',stress));
        xlabel('dissimilarities')
        ylabel('transformed distances')
        axis square
end