%Script to calculate different label ratios from the simulations.

%Luciana de Oliveira, December 2016.

sourceRoot = '/project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/09/092914/targetISruns/';
saveRoot= '/project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/09/092914/targetISruns/';

%Define strings for directory hierarchy as needed
rDDir = {'rD12'};%,'rD60','rD80','rD120','rD140','rD160'};
aPDir = {'aP0p6'};%,,'aP0p7'
outDirNum =1:10;
labelRatio = 0.2;
intensityQuantum=[1 0.3];

%Define number of label ratio

fprintf('\n===============================================================');

%The top level directory is that of receptor density
for  rDDirIndx = 1 : length(rDDir)
    
    tic    
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
      
        for outDirIndx = 1 : length(outDirNum)
            
            %name of current directory
            currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx))];
            
            %Load receptorInfoAll

            tempRecepInfo = load([currDir,filesep,...
                'receptorInfoAll',int2str(outDirNum(outDirIndx)),'.mat']);
            
          receptorInfoAll=tempRecepInfo.receptorInfoAll;
          
          %calculate the receptor info labeled for the different values of
          %label ratio
          
            
          
               receptorInfoLabeled = genReceptorInfoLabeled(receptorInfoAll,...
               labelRatio,intensityQuantum);
          
                
                currOutDir = [saveRoot,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx))];
                
                %Create the direcotry
                mkdir(currOutDir)
                
                %save
             
                save([currOutDir,filesep,'receptorInfoLabeled',int2str(outDirNum(outDirIndx))],'receptorInfoLabeled','-v7.3');
            
                end
                                
                fprintf('... done.');
                
                 end %for each labelRatio
            end
       
        
                            
    %for each aP
    
    elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    


fprintf('\n\nAll done.\n');

clear