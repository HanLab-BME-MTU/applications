function out = bilobePlot(inputData, dataName)
%BILOBEPLOT is a plotting function for the analysis of bilobe patterns
%
% INPUT inputData either n-by-4 array of [spindleLength, bin#, weight,
%                   movie#] from bilobedDistribution, or a cell array with
%                   {spindleLength, projectedIntensities} from
%                   bilobesGfpKinetochore.
%
%



% triage input
if iscell(inputData)
    % we have intensity measurements
    spindleLength = inputData{1};
    intensities = true;
else
    % we have positions
    spindleLength = inputData(:,1);
    intensities = false;
end

if nargin == 1 || isempty(dataName)
    dataName = ' ';
end

plotSigma = true;
plotTrapezoid = true;

% make 26 bins
xLabels = -1/48:1/24:49/48;

% sample every 25 nanometer spindle length
boundaries = [0.7:0.025:1.7;0.9:0.025:1.9];

meanBoundaries = mean(boundaries,1);
nBoundaries = size(boundaries,2);
nBins = length(xLabels);
xall = repmat(xLabels',1,nBoundaries);

% plotData: data needed for plotting
out(1:2) = struct('xall',[],'yall',[],'zall',[],'yTickLabels',[]);

% map everything onto half-spindle
if intensities
    spindleLength = repmat(spindleLength,2,1);
end

% plot averages
figureNames = {'no norm','sum = 1','max = 1'};
for i=1:3

    yall=zeros(nBins,nBoundaries);
    zall=zeros(nBins,nBoundaries);
    zallS = zeros(nBins,nBoundaries);
    yTickLabels = cell(nBoundaries,1);
    nSpindles = zeros(nBoundaries,1);
    spindleIdx = cell(nBoundaries,1);

    % read and potentially symmetrize data
    if intensities
        pInt = inputData{2};
    else
        % dist: bin, weight, movie
        dist = inputData(:,2:end);
    end

    if intensities
        pInt = [pInt(1:13,:),pInt(26:-1:14,:)];
    else
        otherSideIdx = dist(:,1)>13;
        dist(otherSideIdx,1) = 27-dist(otherSideIdx,1);
    end

    % figure, hold on
    for ct = 1:nBoundaries,
        spindleIdx{ct} = (spindleLength>boundaries(1,ct) & spindleLength<boundaries(2,ct));
        if any(spindleIdx{ct})
            if intensities
                % make weighted average - points 0.1um away have no weight
                weights = (0.1-abs(spindleLength(spindleIdx{ct})-meanBoundaries(ct)))/0.1;
                [averageMP,dummy,sigmaMP] = weightedStats(pInt(:,spindleIdx{ct})',...
                    weights,'w');
            else
                % for every bin, sum up the weights multiplied by the
                % distance weights

                % calculate all the weights already now - we need it
                % for the calculatio of n later
                weights = (0.1-abs(spindleLength(spindleIdx{ct})-...
                    meanBoundaries(ct)))/0.1 .*...
                    dist(spindleIdx{ct},2);
                for bin = 13:-1:1
                    % "average": sum the weights in each bin
                    weightsList = weights(dist(spindleIdx{ct},1)==bin);
                     averageMP(bin) = sum(weightsList);
                     sigmaMP(bin) = NaN;
                    % potentially, we can extract distributions for
                    % individual movies, and average those to get a std in
                    % the future, for now, don't have std

                    % two alternatives: (1) Calculate variance between four
                    % areas of equal weight under the triangle filter -
                    % this measures spatial heterogeneity
                    % (2) Bootstrap the variance (take 1000 samples allowing
                    % repetitions) - neat, but what would this really be
                    % measuring?
%                     wDist = bootstrp(1000,'sum',weightsList);
%                     averageMP(bin) = mean(wDist);
%                     sigmaMP(bin) = std(wDist);
                end
                
                

                
            end

            % we will plot (avg,sigma) in two halves of the plot.
            % Therefore remove the other half of both vectors so that we
            % can directly stick them together afterwards

            averageMP = averageMP(1:13);
            sigmaMP = sigmaMP(1:13);

            switch i
                case {3}
                    sigmaMP = sigmaMP/max(averageMP);
                    averageMP = averageMP/max(averageMP);
                case 1
                    % no norming

                case 2
                    sigmaMP = sigmaMP/sum(averageMP);
                    averageMP = averageMP/sum(averageMP);

            end
            %plot3(x,ct*ones(size(x)),z),
            nSpindles(ct) = sum(weights);
            %old: nSpindles(ct) = nnz(spindleIdx{ct});

            zall(:,ct)=[averageMP,averageMP(end:-1:1)];
            zallS(:,ct) = [sigmaMP,sigmaMP(end:-1:1)];
        else
            % there are no spindles this long
            zall(:,ct) = NaN;
            zallS(:,ct) = NaN;
        end
        yall(:,ct)=meanBoundaries(ct);
        yTickLabels{ct}=sprintf('%1.1f/%1.2f', ...
            meanBoundaries(ct),nSpindles(ct));
        
        % convolve with psf if not intensties
        if ~intensities
            dx=xall(2,1)-xall(1,1);
            [FT_XY, FT_Z] = calcFilterParms(0.525,1.4,1.51,'gauss',[], [dx*yall(1,ct) dx*yall(1,ct)]);
            g=gauss1d(-5:5,FT_XY);
            zg=conv(zall(:,ct),g);
            zall(:,ct)=zg(6:end-5);
        end

    end

    % check zall for inf
    zall(~isfinite(zall)) = NaN;
    zallS(~isfinite(zallS)) = NaN;
    % check zall for <0
    zall(zall<0) = 0;

    figure('Name',[dataName,' ',figureNames{i}])
    if plotSigma %&&  ~all(isnan(zallS(:)))
        ah = subplot(1,2,1);
    else
        ah = gca;
    end
    % adjust xall if we want to have a trapezoid plot
    if plotTrapezoid
        xFactor = yall;
        xSub = 0.5;
    else
        xFactor = ones(size(xall));
        xSub = 0;
    end
    contourf((xall-xSub).*xFactor,yall,zall,'LineStyle','none','LevelList',linspace(0,nanmax(zall(:)),100));
    %     figure('Name',figureNames{i}),surf(xall,yall,zall,'FaceColor','interp','FaceLighting','phong')
    %     axis tight
    set(ah,'yTick',meanBoundaries(1:4:end),'yTickLabel',yTickLabels(1:4:end),'yGrid','on','Color',get(gcf,'Color'))

    % set new colormap (see:
    % http://www.research.ibm.com/people/l/lloydt/color/color.HTM;
    % and:
    % http://www.research.ibm.com/dx/proceedings/pravda/index.htm
    % for details)
    cm = isomorphicColormap('b/y');
    colormap(cm)
    
    switch i
        case {3}
            set(ah,'CLim',[0,1])
        case 2
            goodZ = nSpindles>1;
            set(ah,'CLim',[0,nanmax(nanmax(zall(:,goodZ)))])
        case 1
            % here it depends how we normalized before. Implement later, do
            % 01 for now
            if nanmax(zall(:)) < 0.8 || nanmax(zall(:)) > 1
                set(ah,'CLim',[0,nanmax(zall(:))])
            else
                set(ah,'CLim',[0,1])
            end
    end
    if plotSigma %&&  ~all(isnan(zallS(:)))
        % don't put the colorbar here already
    else
        colorbar('peer',ah)
    end
    if plotTrapezoid
        % add lines
        hold on
        for d=0.2:0.2:0.8
            line(d - meanBoundaries(meanBoundaries>=2*d)/2,...
                meanBoundaries(meanBoundaries>=2*d),'Color','k','LineWidth',1);
            line(-d + meanBoundaries(meanBoundaries>=2*d)/2,...
                meanBoundaries(meanBoundaries>=2*d),'Color','k','LineWidth',1);
        end
    else
        % add lines
        hold on
        for d=0.2:0.2:0.8
            line(d./meanBoundaries(meanBoundaries>=2*d),...
                meanBoundaries(meanBoundaries>=2*d),'Color','k','LineWidth',1);
            line(1-d./meanBoundaries(meanBoundaries>=2*d),...
                meanBoundaries(meanBoundaries>=2*d),'Color','k','LineWidth',1);
        end
    end

    if plotSigma %&& ~all(isnan(zallS(:)))
        ah = subplot(1,2,2);
        contourf((xall-xSub).*xFactor,yall,zallS,'LineStyle','none','LevelList',linspace(0,nanmax(zall(:)),100));
        %     figure('Name',figureNames{i}),surf(xall,yall,zall,'FaceColor','interp','FaceLighting','phong')
        %     axis tight
        set(ah,'yTick',meanBoundaries(1:4:end),'yTickLabel',yTickLabels(1:4:end),'yGrid','on','Color',get(gcf,'Color'))

        switch i
            case {3}
                set(ah,'CLim',[0,1])
            case 2
                goodZ = nSpindles>1;
                set(ah,'CLim',[0,nanmax(nanmax(zall(:,goodZ)))])
            case 1
                % here it depends how we normalized before. Implement later, do
                % 01 for now
                if nanmax(zall(:)) < 0.8 || nanmax(zall(:)) > 1
                    set(ah,'CLim',[0,nanmax(zall(:))])
                else
                    set(ah,'CLim',[0,1])
                end
        end
        colorbar('peer',ah)
        % add lines
        if plotTrapezoid
            % add lines
            hold on
            for d=0.2:0.2:0.8
                line(d - meanBoundaries(meanBoundaries>=2*d)/2,...
                    meanBoundaries(meanBoundaries>=2*d),'Color','k','LineWidth',1);
                line(-d + meanBoundaries(meanBoundaries>=2*d)/2,...
                    meanBoundaries(meanBoundaries>=2*d),'Color','k','LineWidth',1);
            end
        else
            hold on
            for d=0.2:0.2:0.8
                line(d./meanBoundaries(meanBoundaries>=2*d),...
                    meanBoundaries(meanBoundaries>=2*d),'Color','k','LineWidth',1);
                line(1-d./meanBoundaries(meanBoundaries>=2*d),...
                    meanBoundaries(meanBoundaries>=2*d),'Color','k','LineWidth',1);
            end
        end
    end


    out(i).xall = xall;
    out(i).yall = yall;
    out(i).zall = zall;
    out(i).zallS = zallS;
    out(i).yTickLabels = yTickLabels;
    out(i).spindleIdx = spindleIdx;
    out(i).nSpindles = nSpindles;


end

% plot distributions
figure('Name','Distributions')
subplot(2,2,1)
hold on
cmap = jet(nBoundaries);
for i=1:nBoundaries
    plot(out(1).xall(1:13,1).*out(1).yall(1:13,i),out(1).zall(1:13,i),'Color',cmap(i,:));
end
for i=1:3
subplot(2,2,i+1)
plot(out(i).yall(13,:),out(i).zall(13,:))
end



% if intensities
%     % adapt for positions later
%     figure('Name',sprintf('%s individual; sum=1',dataName));
%     [sortedSpindleLength,sortIdx]=sort(inputData{1});
%     int = inputData{2};
%     int = int(:,sortIdx)';
%     int = int./repmat(sum(int,2),1,size(int,2));
%     intSymm = 0.5*(int + int(:,end:-1:1));
%     subplot(1,2,1),imshow([int,NaN(size(int,1),1),intSymm],[]),
%     colormap jet,
%     subplot(1,2,2),
%     plot(sortedSpindleLength,length(sortIdx):-1:1,'-+')
% end
%
% % loop to make histograms
% stages = [1,1.2,1.6,2];
%
% % read data. symmetrize distances
% if intensities
%     pInt = inputData{2};
% else
%     % dist: bin, weight, movie
%     dist = inputData(:,2:end);
%     % reverse bins
%     dist(2:2:end,1) = 27-dist(2:2:end,1);
% end
% % allMP = (allMP+allMP(end:-1:1,:))/2;
%
% for ct = 1:length(stages)-1,
%     figure('Name',sprintf('%s %1.1f -> %1.1f',dataName, stages(ct:ct+1)));
%     sidx = find(spindleLength>stages(ct) & spindleLength<stages(ct+1));
%     if intensities
%         averageMP = mean(pInt(:,sidx),2);
%     else
%         % for every bin, sum up the weights
%
%         % calculate all the weights already now - we need it
%         % for the calculatio of n later
%         weights = dist(spindleIdx{ct},2);
%         for bin = 26:-1:1
%             % "average": sum the weights in each bin
%             averageMP(i) = sum(weights(dist(spindleIdx{ct},1)==bin));
%         end
%     end
%
%     samp=sum(averageMP);
%     averageMP = averageMP/samp;
%
%     hold on
%     if intensities
%         for i=1:length(sidx)
%             plot(xLabels,pInt(:,sidx(i))/samp,'b');
%         end
%     else
%         % do individual movies later
%     end
%     plot(xLabels,averageMP,'r','LineWidth',1.5)
% end