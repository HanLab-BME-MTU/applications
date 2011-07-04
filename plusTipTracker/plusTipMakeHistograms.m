function plusTipMakeHistograms(speedLifeDispMat,saveDir)
% plusTipMakeHistograms saves speed, lifetime, and displacement histograms
% for growth, fgap, and bgap populations
%
% Called by plusTipPoolGroupData.  See function for input format.
% Kathryn Applegate, Jan 2010
% Sebastien Besson, June 2011

if ~isdir(saveDir), mkdir(saveDir); end

% Define types of data and events
dataType(1).name = 'speed';
dataType(1).xlabel = 'Speed (microns/min)';
dataType(1).columns = 1:3;
dataType(2).name = 'lifetime';
dataType(2).xlabel = 'Lifetime (sec)';
dataType(2).columns = 4:6;
dataType(3).name = 'displacement';
dataType(3).xlabel = 'Displacement (microns)';
dataType(3).columns = 7:9;

eventType(1).name='growth';
eventType(1).color=[1 0 0];
eventType(2).name='fgap';
eventType(2).color=[0 0 1];
eventType(3).name='bgap';
eventType(3).color=[0 1 0];

for i=1:numel(dataType)
    data = speedLifeDispMat(:,dataType(i).columns);

    % create x-axis bins spanning all values
    n=linspace(nanmin(data(:)),nanmax(data(:)),25);

    % bin the samples
    binData = arrayfun(@(x)histc(data(:,x),n),1:3,'UniformOutput',false);

    % put the binned values into a matrix filled with NaNs for the stacked plot 
    maxLength = max(cellfun(@length,binData));
    M = cell2mat(cellfun(@(x) vertcat(x,NaN(maxLength-length(x),1)),...
        binData,'UniformOutput',false));

    % Make the stacked plot
    figure
    bar(n,M,'stack')
    colormap(vertcat(eventType.color))
    legend({eventType.name},'Location','best')
    title(['Stacked ' dataType(i).name ' distributions'])
    xlabel(dataType(i).xlabel);
    ylabel('Frequency of tracks');
    saveas(gcf,[saveDir filesep 'histogram_' dataType(i).name '_stacked.tif'])
    close(gcf)

    % Make individual plots for non-empty event data
    validBinData = find(~cellfun(@isempty,binData));
    for j=validBinData
        figure;
        bar(n,binData{j})
        colormap(eventType(j).color)
        title([eventType(j).name ' ' dataType(i).name ' distribution'])
        xlabel(dataType(i).xlabel);
        ylabel('Frequency of tracks');
        saveas(gcf,[saveDir filesep 'histogram_' dataType(i).name '_'...
            eventType(j).name '.tif'])
        close(gcf);
    end
end