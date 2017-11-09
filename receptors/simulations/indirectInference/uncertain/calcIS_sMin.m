function calcIS_sMin(pathToS)
%CALCIS_SMIN collects the minimum S values (Mahalanobis distances) from the
%set of s values in sMatrix.  This function is part of the indirect 
%inference based model calibration framework.
%
%   INPUT:
%           pathToS:    location of sMatrix mat files
%
%   OUTPUT:
%           the set of s minima are saved in a mat file called sMins. The
%           minima are saved in a matrix with the following columns:
%           col. 1:     probe receptor densities
%           col. 2:     association prob with minimum S for the density in
%                       col. 1
%           col. 3:     label ratio with minimum S for the density in
%                       col. 1
%           col. 4:     the minimum S value
%
%   Robel Yirdaw, 11/25/14
%

    %These are currently the set of densities, association probablities and
    %label ratios for probe intermediate statistics
    rDStr = {'rD2';'rD4';'rD6';'rD8';'rD10';'rD12';'rD14';'rD16'};
    rDvals = [2;4;6;8;10;12;14;16];
    aPvals = [0.2;0.3;0.4;0.5;0.6;0.7;0.8];
    lRvals = [0.1;0.2;0.3;0.4;0.5;0.6];
    
    %An array with three columns to hold the minimum values for each
    %receptor density set. Columns: receptor density, assocProb, 
    %labelRatio and the minimum S value
    
    minSvals = NaN(length(rDStr),4);
    minSvals(:,1) = rDvals;
    
    for rDindx=1:length(rDStr)
        
        loadStruct = load([pathToS,'sMatrix_',rDStr{rDindx},'.mat']);
        
        %Determine the minimum S value and populate columns at current
        %density row
        minSvals(rDindx,4) = min(min(loadStruct.sMatrix));
        [row,col,~] = find(loadStruct.sMatrix == minSvals(rDindx,4));
        minSvals(rDindx,2:3) = [aPvals(col) lRvals(row)];
                
        clear loadStruct row col
        
    end % for each receptor density
    
    %Write to file
    outFile = [pathToS,'sMins.mat'];
    save(outFile,'minSvals');
    
    fprintf('\nMinimum S values saved in %s.\n',outFile);
    
end % function