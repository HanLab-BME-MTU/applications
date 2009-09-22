function plusTipMakeHistograms(speedLifeDispMat,saveDir)
% plusTipMakeHistograms saves speed, lifetime, and displacement histograms
% for growth, fgap, and bgap populations

if ~isdir(saveDir)
    mkdir(saveDir)
end


for iParam=1:3
    switch iParam
        case 1
            data=speedLifeDispMat(:,1:3); % speed
            titleStr='speed';
            xStr='speed (microns/min)';
        case 2
            data=speedLifeDispMat(:,4:6); % lifetime
            titleStr='lifetime';
            xStr='lifetime (min)';
        case 3
            data=speedLifeDispMat(:,7:9); % displacement
            titleStr='displacement';
            xStr='displacement (microns)';
    end


    % create x-axis bins spanning all values
    n=linspace(nanmin(data(:)),nanmax(data(:)),25);

    % bin the samples
    [x1,dummy] = histc(data(:,1),n); % growth
    [x2,dummy] = histc(data(:,2),n); % fgap
    [x3,dummy] = histc(data(:,3),n); % bgap

    % put the binned values into a matrix for the stacked plot
    M=nan(max([length(x1) length(x2) length(x3)]),3);
    M(1:length(x1),1)=x1;
    M(1:length(x2),2)=x2;
    M(1:length(x3),3)=x3;

    % make the plot
    figure
    bar(n,M,'stack')
    colormap([1 0 0; 0 0 1; 0 1 0])
    legend('growth','fgap','bgap','Location','best')
    title(['stacked ' titleStr ' distributions'])
    xlabel(xStr);
    ylabel('frequency of tracks');
    %saveas(gcf,[saveDir filesep titleStr '_stackedHist.fig'])
    saveas(gcf,[saveDir filesep titleStr '_stackedHist.tif'])
    close(gcf)

    figure;
    % growth
    if ~isempty(x1)
        bar(n,x1,'r')
        title(['growth ' titleStr ' distribution'])
        xlabel(xStr);
        ylabel('frequency of tracks');

        %saveas(gcf,[saveDir filesep titleStr '_growthHist.fig'])
        saveas(gcf,[saveDir filesep titleStr '_growthHist.tif'])
    end
    close(gcf)

    figure
    % fgap
    if ~isempty(x2)
        bar(n,x2,'b')
        title(['fgap ' titleStr ' distribution'])
        xlabel(xStr);
        ylabel('frequency of tracks');

        %saveas(gcf,[saveDir filesep titleStr '_fgapHist.fig'])
        saveas(gcf,[saveDir filesep titleStr '_fgapHist.tif'])
    end
    close(gcf)

    figure
    % bgap
    if ~isempty(x3)
        bar(n,x3,'g')
        title(['bgap ' titleStr ' distribution'])
        xlabel(xStr);
        ylabel('frequency of tracks');

        %saveas(gcf,[saveDir filesep titleStr '_bgapHist.fig'])
        saveas(gcf,[saveDir filesep titleStr '_bgapHist.tif'])
    end
    close(gcf)

end