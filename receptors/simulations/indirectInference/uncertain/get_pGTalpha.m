function [pGTalpha] = get_pGTalpha(pathToP,alpha)
%GET_PGTALPHA collects p-values that are greater than alpha.  This function
%is part of the indirect inference based model calibration framework.
%
%   INPUT:
%           pathToP:    location of p-values
%           alpha:      threshold value
%
%   OUTPUT:
%           pGTalpha:   p-values that are greater than alpha arranged in a
%                       cell for each probe receptor density.  Each cell
%                       element contains a matrix with four columns when a
%                       p-value greater than alpha is found, ([] otherwise)
%                       col 1:  the receptor density
%                       col 2:  assoc. prob for the corresponding p-value
%                       col 3:  label ratio for the corresponding p-value
%                       col 4:  the p-value that is > alpha
%
%   Robel Yirdaw, 11/25/14
%

    %These are currently the set of densities, association probablities and
    %label ratios for probe intermediate statistics
    rDStr = {'rD2';'rD4';'rD6';'rD8';'rD10';'rD12';'rD14';'rD16'};
    rDvals = [2;4;6;8;10;12;14;16];
    aPvals = [0.2;0.3;0.4;0.5;0.6;0.7;0.8];
    lRvals = [0.1;0.2;0.3;0.4;0.5;0.6];
    
    %Initialize the cell vector for p-values and related info to be
    %collected
    pGTalpha = cell(length(rDStr),1);
    
    %Iterate through each density since pMatrix is saved by density
    for rDindx=1:length(rDStr)
        
        %Load pMatrix
        loadStruct = load([pathToP,'pMatrix_',rDStr{rDindx},'.mat']);
        tempPmatrix = loadStruct.pMatrix;
        %Determine the p-values that are > alpah and populate columns 
        %at current density row        
        [row,col,~] = find(tempPmatrix > alpha);
        pGTalpha{rDindx} = [repmat(rDvals(rDindx),length(row),1) aPvals(col)...
            lRvals(row) tempPmatrix(sub2ind(size(tempPmatrix),row,col))];
        
        clear loadStruct tempPmatrix row col
    end
    
    %Write to file
    outFile = [pathToP,'pGTalpha.mat'];
    save(outFile,'pGTalpha');
    
    fprintf('\nP-values > %d saved in %s.\n',alpha,outFile);
    
end
