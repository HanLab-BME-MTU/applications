function getMultipleAggregStateLabeldReceptorsBoot(sourceRoot,saveRoot,rDDir,dRDir,aPDir,outDirNum,lRDir,labelRatio,bootRepSize,bootsCond)
%function to calculate aggregation state from receptorInfoAll, for 100% of receptors, for the static
%data
% SYNOPSIS functionCalcMultipleAggregStateFromReceptorInfoAlldR(sourceRoot,saveRoot,rDDir,aPDir,dRDir,outDirNum,lRDir)
%
%INPUT
%  sourceRoot : path for the receptorInfoLabeled.
%
%  saveRoot   : path where the detection will be saved.
%
%  rDDir: multiple receptor density directories that will be analyzed. It
%  is a cell in the form: {'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'}
%
%  aPDir      : association probability directories. A cell in the form: {'aP0p2','aP0p3','aP0p4'}
%
%  dRDir      : dissociation rate directories. A cell in the form: {'aP0p2','aP0p3','aP0p4'}
%
%  outDirNum  : list with the simulation number. 1:10 for the low density and
%  1:30 for the high density.
%
%  lRDir      : labed receptor directory. A cell in the form: {'lR0p2';'lR0p3';'lR0p4'}
%         
%  labeled ratio: fraction of labeled receptors.
%  bootRepSize: number of repetitions boot
%
%  bootsCond: 0 for no bootstrapping, 1 for bootstrapping
%
% OUTPUT
% There is no output, the results will be saved in the sourceRoot directory
% as detectionAggregState.m
%
% Luciana de Oliveira, May 2017

fprintf('\n===============================================================');

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    for dRDirIndx = 1 : length(dRDir)
        tic
        %Iterate through association probability values per density
        for aPDirIndx = 1 : length(aPDir)
            
            fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
            
            %iterate through the different runs
            for outDirIndx = 1 : length(outDirNum)
                
                for lRDirIndx = 1 : length(lRDir)
                    
                    %name of current directory
%                    currDir = [originRoot,filesep,rDDir{rDDirIndx},filesep, aPDir{aPDirIndx},filesep,...
%                        'out',int2str(outDirNum(outDirIndx))];
   currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,  dRDir{dRDirIndx},filesep, aPDir{aPDirIndx},filesep,...
                              'out',int2str(outDirNum(outDirIndx))];
%                     
                    %load compTracks
                    tmp = load([currDir,filesep,'receptorInfoAll' int2str(outDirNum(outDirIndx)) '.mat']);
                    receptorInfoAll = tmp.receptorInfoAll;
                    clear tmp
                    %reserve space for the detectionAggregState for all the
                    %repetitions of bootstrapping
                    
                    detectionAggregState=cell(bootRepSize,1);
                    
                    
                    for iBoot=1:bootRepSize
                    
                                        %get aggregation state
                    detectionAggregState {iBoot}= genAggregStateFromReceptorLabeledSuperResBoots(receptorInfoAll,labelRatio,bootsCond);
                    end
                    fprintf('\n   Out = %d',outDirIndx);
                    
                    %save tracks with aggregation state
                    saveDir = [saveRoot,filesep,rDDir{rDDirIndx},filesep,dRDir{dRDirIndx},...
                        filesep,aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx)),...
                        filesep,lRDir{lRDirIndx}];
                    
                mkdir(saveDir)
                    
                    save([saveDir,filesep,'detectionAggregState'],'detectionAggregState','-v7.3');
                    clear detectionAggregState
                    
                    fprintf('... done.');
                end
                
            end %for each outDir
            
        end %for each aP
    end
    elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    
end %for each rD

fprintf('\n\nAll done.\n');

clear
