%script for check if the simulations in the same condition are following a
%trend. It is done after the determination of cluster population for
%superresolution.


sourceRoot='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170424/probeISruns';
saveDir='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170502/simulationswithproblems';

%Define strings for directory hierarchy as needed
rDDir = {'rD12'}; %{'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'};%,'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'
numberReceptorsInput=12;
aPDir = {'aP0p4','aP0p6'}; %{'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};%'aP0p6','aP0p7','aP0p8' ,'aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'
dRDir={'dR0p5','dR1p5'};%,'dR1p5','dR1p75'
outDirNum =1:10;
lRDir = {'lR1'};
intensityQuantum=[1 0.3];


outProblem=cell(length(rDDir),length(aPDir),length(dRDir));
for rDDirIndx = 1 : length(rDDir)
    for dRDirIndx = 1 : length(dRDir)
        tic
        %Iterate through association probability values per density
        for aPDirIndx = 1 : length(aPDir)
            
            fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
            
            %iterate through the different labeling ratios
            for lRDirIndx = 1 : length(lRDir)
                
                %name of current directory
                
                %define output cell array
                densityPerMovie = zeros(5,length(outDirNum));
                
                %iterate through the different runs
                for outDirIndx = 1 : length(outDirNum)
%                      currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep, aPDir{aPDirIndx},filesep,  dRDir{dRDirIndx},filesep,...
%                               'out',int2str(outDirNum(outDirIndx))];
                    
%                    currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,  dRDir{dRDirIndx},filesep, aPDir{aPDirIndx},filesep,...
%                               'out',int2str(outDirNum(outDirIndx))];
                    
                          currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep, aPDir{aPDirIndx},filesep,...
                              'out',int2str(outDirNum(outDirIndx))];
                    
                          
                    tmp = load([currDir,filesep,'receptorInfoAll' int2str(outDirNum(outDirIndx)) '.mat']);
                    receptorInfoAll = tmp.receptorInfoAll;
                    clear tmp
                   
                                       %store in cell array
                  receptorPerCluster = receptorInfoAll.clust2receptAssign(:,:,end); 
                   % calculate the number of receptors per cluster size
maxClusterSize=size(receptorPerCluster,2);
                    % reserve space for the detectionAggregState
numberOfReceptors=zeros(1,maxClusterSize);

% calculate the density of receptors in each cluster size

for i=1:maxClusterSize
    numberOfReceptors(i) = nnz(receptorPerCluster(:,i));
end
%check if the number of receptors is correct
numberReceptorSim=sum(numberOfReceptors)/areaSideLen^2;
if numberReceptorSim~=numberReceptorsInput(rDDirIndx)
  outProblem{rDDirIndx,aPDirIndx,dRDirIndx}=currDir;  
%    mkdir(saveDir)
%                     save([saveDir,'/outProblem'],'outProblem','-v7.3'); 
end
                end %for each outDir
                              % if there is a large difference between these values, save the
                % path
                
              
                
            end %for each labelRatio
            
        end %for each aP
        
        
    end
    
    elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    
end %for each rD


fprintf('\n\nAll done.\n');