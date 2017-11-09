%Script to calculate aggregation state from receptorInfoAll, for 100% of receptors, for the static
%data
%
%Luciana de Oliveira, February 2017

sourceRoot = '/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170317/varyDissRate/probe';
saveRoot='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170516/superRes';
% saveRoot='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170403/superRes/target';


%Define strings for directory hierarchy as needed
rDDir ={'rD4'};%,'rD20','rD40','rD60','rD80','rD100','rD120','rD140''rD8',
aPDir = {'aP0p5'};%,'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'}
outDirNum =1:10;
dRdir={'dR0p75'};
lRDir = {'lR1'};%{'lR0p14';''lR0p1';'lR0p2';'lR0p3';'lR0p4';'lR0p5';'lR0p6'

fprintf('\n===============================================================');

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
    tic
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
        
        fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
        
        %iterate through the different runs
        for outDirIndx = 1 : length(outDirNum)
            
              for lRDirIndx = 1 : length(lRDir)
            
            %name of current directory
            currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx))];
            
            %load compTracks
            tmp = load([currDir,filesep,'receptorInfoAll' int2str(outDirNum(outDirIndx)) '.mat']);
            receptorInfoAll = tmp.receptorInfoAll;
            clear tmp
            
            %get aggregation state
            detectionAggregState = aggregStateFromReceptorInfoAll(receptorInfoAll);
            
            fprintf('\n   Out = %d',outDirIndx);
            
            %save tracks with aggregation state
            saveDir = [saveRoot,filesep,rDDir{rDDirIndx},filesep,dRdir{1},filesep, ...
                aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx)),...
                    filesep,lRDir{lRDirIndx}];
                   mkdir(saveDir)
            
            save([saveDir,filesep,'detectionAggregState'],'detectionAggregState','-v7.3');
            clear detectionAggregState
            
            fprintf('... done.');
              end
              
        end %for each outDir
        
    end %for each aP
    
    elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    
end %for each rD

fprintf('\n\nAll done.\n');

clear
