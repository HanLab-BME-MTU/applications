%Script to calculate aggregation state from receptorInfoAll, for 100% of receptors, for the static
%data
%
%Luciana de Oliveira, February 2017

sourceRoot = '/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170426/varyDissRate/probeISruns';
saveRoot='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170425/varyDissRate/probe/probeIS_sT25_dT0p1';

%Define strings for directory hierarchy as needed
 rDDir = {'rD6'};%,'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'
 aPDir =  {'aP0p5'};%'aP0p6','aP0p7','aP0p8' ,'aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'
 dRDir={'dR0p25'};%,'dR1p5','dR1p75'
outDirNum =1:10;
lRDir = {'lR1'};%{'lR0p14';''lR0p1';'lR0p2';'lR0p3';'lR0p4';'lR0p5';'lR0p6'

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
                   currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,dRDir{dRDirIndx},filesep, aPDir{aPDirIndx},filesep,...
                       'out',int2str(outDirNum(outDirIndx))];
%   currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep, aPDir{aPDirIndx},filesep,...
%                       dRDir{dRDirIndx},filesep,  'out',int2str(outDirNum(outDirIndx))];
%                     
                    %load compTracks
                    tmp = load([currDir,filesep,'receptorInfoAll' int2str(outDirNum(outDirIndx)) '.mat']);
                    receptorInfoAll = tmp.receptorInfoAll;
                    clear tmp
                    
                    %get aggregation state
                    detectionAggregState = aggregStateFromReceptorInfoAll(receptorInfoAll);
                    
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
