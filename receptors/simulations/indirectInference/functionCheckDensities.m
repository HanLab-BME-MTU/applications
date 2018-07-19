function simWithDensityProblems=functionCheckDensities(sourceRoot,rDDir,aPDir,dRDir,outDirNum,rDVal,areaSideLen)

%INPUT
%  sourceRoot : path for the out.mat.
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
%  areaSideLen:      Simulation/image side length values,
%                   which can be a single value or a value per
%                   side. In units of interest (e.g. um).
%
% OUTPUT
% There is no output, the results will be saved in the sourceRoot directory
% as compTracks.m
%
% Luciana de Oliveira, May 2017.


lRDir = {'lR1'};
numberReceptorsInput=rDVal;
outProblem=cell(length(rDDir),length(aPDir),length(dRDir));
for rDDirIndx = 1 : length(rDDir)
    for dRDirIndx = 1 : length(dRDir)
        
        %Iterate through association probability values per density
        for aPDirIndx = 1 : length(aPDir)
            
                      
            %iterate through the different labeling ratios
            for lRDirIndx = 1 : length(lRDir)
                
                %name of current directory
                
                %define output cell array
                densityPerMovie = zeros(5,length(outDirNum));
                
                %iterate through the different runs
                for outDirIndx = 1 : length(outDirNum)
                    %                      currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep, aPDir{aPDirIndx},filesep,  dRDir{dRDirIndx},filesep,...
                    %                               'out',int2str(outDirNum(outDirIndx))];
                     fprintf('\nProcessing rD = %s, aP = %s, dR = %s, out=%s ',rDDir{rDDirIndx},aPDir{aPDirIndx}, dRDir{dRDirIndx}, int2str(outDirNum(outDirIndx)));
                     
                                        currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,  dRDir{dRDirIndx},filesep, aPDir{aPDirIndx},filesep,...
                                                 'out',int2str(outDirNum(outDirIndx))];
                    
%                     currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep, aPDir{aPDirIndx},filesep,...
%                         'out',int2str(outDirNum(outDirIndx))];
                    
                    
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
                    end
                   
                end %for each outDir
                
            end %for each labelRatio
            
        end %for each aP
        
        
    end
    
    
end %for each rD

simWithDensityProblems=outProblem;

