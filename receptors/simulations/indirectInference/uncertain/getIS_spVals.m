function [sMatrix,pMatrix] = getIS_spVals(pathToTargetIS)
%GETIS_SPVALS locates intermediate statistics for target and probes,
%invokes calcIS_sVal_cutOff and calcIS_pVal_cutOff to calculate S and
%p-values which are then organized in matrices and output to file.  This 
%function is part of the indirect inference based model calibration 
%framework.
%
%   INPUT:
%           pathToTargetIS:     location of target intermediate statistics
%
%   OUTPUT:
%           sMatrix:            calculated Mahalanobis distance S in a
%                               matrix for each *probe* density and with
%                               columns corresponding to association
%                               probabilities and rows corresponding to
%                               label ratios. Also written to file.
%           pMatrix:            calculated p-values in a
%                               matrix for each *probe* density and with
%                               columns corresponding to association
%                               probabilities and rows corresponding to
%                               label ratios. Also written to file.
%           isVals_mod:         intermediate statistics vector and matrices
%                               may have been modified when calculating S.
%                               The modified quantities are returned by
%                               calcIS_sVal_cutOff as a struct. Save this 
%                               to file.           
%   Robel Yirdaw, 10/02/14
%       Modified, 11/25/14 (probe rD2 and lR0.6)
%

    targetISroot = '/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20160721_Analysis/target/';
    %Location of probe IS other than the two below
    pathToProbeIS1 = '/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20160721_Analysis/probe/';
    %Location of probe lR = 0.6 for rD >= 4
    pathToProbeIS2 = '/project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/10/103014/probeISanalysis_sT25_dT0p01/';
    %Location of probe rD >= 2
    pathToProbeIS3 = '/project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/11/112414/probeISanalysis_sT25_dT0p01/';
    rDDir = {'rD2';'rD4';'rD6';'rD8';'rD10';'rD12';'rD14';'rD16'};
    aPDir = {'aP0p2';'aP0p3';'aP0p4';'aP0p5';'aP0p6';'aP0p7';'aP0p8'};
    lRDir = {'lR0p1';'lR0p2';'lR0p3';'lR0p4';'lR0p5';'lR0p6'};   
    
    sMatrix = NaN(length(lRDir),length(aPDir));
    pMatrix = NaN(length(lRDir),length(aPDir));
    %Will be saving modified IS values.
    isVals_mod = cell(length(lRDir),length(aPDir));    
    
    %Degrees of freedom for p-value calucation
    DOF = 9;    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Load target intermediate statistics values
    %pathToTargetIS = targetDirFile.pathToTargetIS;
    targetIS_fileName = 'isVals_target_cutOff';
    
    loadStruct_target = load([targetISroot,pathToTargetIS,targetIS_fileName]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Echo paths
    fprintf('\n======================================================================================');
    fprintf('\nPath to target IS is %s.',pathToTargetIS);
    fprintf('\nPath to probe IS 1 is %s.',pathToProbeIS1);
    fprintf('\nPath to probe IS 2 is %s.',pathToProbeIS2);
    fprintf('\nPath to probe IS 3 is %s.',pathToProbeIS3);
    fprintf('\n======================================================================================');
    fprintf('\n');
    
    %Iterate through probe receptor density, association probability and
    %label ratio
    for rDDirIndx=1:length(rDDir)
        
        for aPDirIndx=1:length(aPDir)

            for lRDirIndx=1:length(lRDir)

                fprintf('\nProcessing probe IS for %s, %s, %s...',rDDir{rDDirIndx},...
                    aPDir{aPDirIndx},lRDir{lRDirIndx});

                %Load probe intermediate statistics values
                %Probe IS values are now at multiple locations!
                if (rDDirIndx == 1)
                    pathToProbeIS = pathToProbeIS3;
                elseif ((rDDirIndx > 1) && (lRDirIndx == 6))
                    pathToProbeIS = pathToProbeIS2;
                else
                    pathToProbeIS = pathToProbeIS1;
                end
                
                loadStruct_probe = load([pathToProbeIS,rDDir{rDDirIndx},'/',...
                    aPDir{aPDirIndx},'/',lRDir{lRDirIndx},'/','isVals_',...
                    rDDir{rDDirIndx},'_',aPDir{aPDirIndx},'_',lRDir{lRDirIndx},'_cutOff.mat']);

                %Calculate S
                [sMatrix(lRDirIndx,aPDirIndx),isVals_mod{lRDirIndx,aPDirIndx}] =...
                    calcIS_sVal_cutOff(loadStruct_target.isVals,loadStruct_probe.isVals);
                %Calculate p-value
                pMatrix(lRDirIndx,aPDirIndx) = calcIS_pVal_cutOff(sMatrix(lRDirIndx,aPDirIndx),DOF);            

                fprintf('done.');
            end % for each labelRatio

        end % for each assocProb
        
        %Set up output files
        isVals_mod_outFile = [targetISroot,'/','isVals_mod_',rDDir{rDDirIndx},'.mat'];                
        sMatrixOutFile = [targetISroot,'sMatrix_',rDDir{rDDirIndx},'.mat'];
        pMatrixOutFile = [targetISroot,'pMatrix_',rDDir{rDDirIndx},'.mat'];

        %Write values to files
        save(isVals_mod_outFile,'isVals_mod');        
        save(sMatrixOutFile,'sMatrix');
        save(pMatrixOutFile,'pMatrix');

        fprintf('\nModified IS values saved in %s.',isVals_mod_outFile);                
        fprintf('\ns matrix saved in %s.',sMatrixOutFile);
        fprintf('\np matrix saved in %s.',pMatrixOutFile);
        fprintf('\n\nDone.');
        
    end % for each receptor density
    
end %function
    
    
            
            
    