function plotIS_pMatrix(pathToP)
%PLOTIS_PMATRIX plot p-values saved in pMatrix. This function is part of 
%the indirect inference based model calibration framework.
%
%   INPUT:
%           pathToP:        loation of pMatrix
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
    
    %Iterate through each set of p-values, saved for each probe density
    for rDindx=1:length(rDStr)
        %Load matrix
        loadStruct = load([pathToP,'pMatrix_',rDStr{rDindx},'.mat']);
        %Plot p-values
        imagesc(aPvals,lRvals,loadStruct.pMatrix);
        %Adjust axes and labels
        figH = gcf;
        set(figH,'Name',rDStr{rDindx});
        
        axH = gca;
        set(axH,'YDir','normal');
        set(axH,'FontSize',12);
        xlabel(axH,'Association probability','FontSize',12);        
        ylabel(axH,'Label ratio','FontSize',12);
        title(axH,['Target ',pathToP,', probe ',rDStr{rDindx}],'FontSize',12);
            
        %Save figure
        outFile = [pathToP,'pMatrix_',rDStr{rDindx},'_plot'];
        saveas(figH,outFile,'png');
        saveas(figH,outFile,'fig');
        
        fprintf('\nFigures saved in %s.\n',outFile);
        
    end % for each receptor density      
    
    
end % function

