function plotIS_sMatrix(pathToS)
%PLOTIS_SMATRIX plots S values saved in sMatrix. This function is part of 
%the indirect inference based model calibration framework.
%
%   INPUT:
%           pathToS:        loation of sMatrix
%
%   OUTPUT:
%           none.
%
%   Robel Yirdaw, 11/25/14
%

    %These are currently the set of densities, association probablities and
    %label ratios for probe intermediate statistics
    rDStr = {'rD2';'rD4';'rD6';'rD8';'rD10';'rD12';'rD14';'rD16'};
    %rDvals = [4;6;8;10;12;14;16];
    aPvals = [0.2;0.3;0.4;0.5;0.6;0.7;0.8];
    lRvals = [0.1;0.2;0.3;0.4;0.5;0.6];
    
    %Iterate through each set of S values, saved for each probe density    
    for rDindx=1:length(rDStr)
        
        loadStruct = load([pathToS,'sMatrix_',rDStr{rDindx},'.mat']);
        %Create figure and set its NextPlot to 'add'
        figH = figure('Name',rDStr{rDindx});
        axH = gca;
        set(axH,'NextPlot','add');
        
        %Determine a normalizing factor for bubble sizes, based on a bubble
        %of size 5000
        normFactor = max(max(loadStruct.sMatrix))/5000;
        
        %Normalize values
        tempS = loadStruct.sMatrix/normFactor;
        
        %Begin plotting
        for lRindx=1:length(lRvals)
            scatter(aPvals,repmat(lRvals(lRindx),length(aPvals),1),tempS(lRindx,:));
        end
        
        %Adjust axes and labels
        xlim(axH,[0 1]);
        ylim(axH,[0 0.6]);

        set(axH,'FontSize',12);
        xlabel(axH,'Association probability','FontSize',12);        
        ylabel(axH,'Label ratio','FontSize',12);
        title(axH,['Target ',pathToS,', probe ',rDStr{rDindx}],'FontSize',12);
            
        %Save figure
        outFile = [pathToS,'sMatrix_',rDStr{rDindx},'_plot'];
        saveas(figH,outFile,'png');
        saveas(figH,outFile,'fig');
        
        fprintf('\nFigures saved in %s.\n',outFile);
        
    end % for each receptor density      
    
    
end % function

