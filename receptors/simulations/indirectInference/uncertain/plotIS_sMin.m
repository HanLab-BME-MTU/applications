function plotIS_sMin(pathToS)
%PLOTIS_SMIN plots minimum S values saved in sMins.mat. This 
%function is part of the indirect inference based model calibration 
%framework.
%
%   INPUT:
%           pathToS:        loation of sMins
%
%   OUTPUT:
%           none.
%
%   Robel Yirdaw, 11/25/14
%

    %These are currently the set of densities
    rDStr = {'rD2';'rD4';'rD6';'rD8';'rD10';'rD12';'rD14';'rD16'};
    
    loadStruct = load([pathToS,'sMins.mat']);
    
    %Create figure and set its NextPlot to 'add'
    figH = figure('Name',pathToS);
    axH = gca;
    set(axH,'NextPlot','add');
        
    %Determine a scaling factor for bubble sizes, based on a bubble
    %of size 5000
    scaleFactor = 5000/max(max(loadStruct.minSvals(:,4)));
        
    %Scale values
    tempMins = loadStruct.minSvals(:,4)*scaleFactor;
        
    %Begin plotting
    for indx=1:length(loadStruct.minSvals(:,1))
        tempH = scatter(loadStruct.minSvals(indx,2),loadStruct.minSvals(indx,3),...,
            tempMins(indx,1));
        set(tempH,'DisplayName',rDStr{indx});        
    end

    %Adjust axes and lables
    xlim(axH,[0 1]);
    ylim(axH,[0 0.7]);

    set(axH,'FontSize',12);
    xlabel(axH,'Association probability','FontSize',12);        
    ylabel(axH,'Label ratio','FontSize',12);
    title(axH,pathToS,'FontSize',12);

    legend(axH,'show');
    legend('boxoff');
    
    %Save figure
    outFile = [pathToS,'sMins_plot'];
    saveas(figH,outFile,'png');
    saveas(figH,outFile,'fig');

    fprintf('\nFigures saved in %s.\n',outFile);
           
end % function

