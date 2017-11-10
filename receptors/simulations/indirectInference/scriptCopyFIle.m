%script to copy files from a directory to another

sourceRoot ='/home2/s169185/Documents/analysisCOpy/all other/analysis/probe';
destinationRoot ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170501';
%Define strings for directory hierarchy as needed




rDDir = {'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'};%,'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'
aPDir =  {'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};%'aP0p6','aP0p7','aP0p8' ,'aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'
dRDir={'dR0p75','dR1','dR1p25','dR1p5','dR1p75'};%,'dR1p5','dR1p75'
outDirNum =1:10;
% lRDir = {'lR1'};%{'lR0p14';''lR0p1';'lR0p2';'lR0p3';'lR0p4';'lR0p5';'lR0p6'

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
                
%                 for lRDirIndx = 1 : length(lRDir)
                    
                    %name of current directory
                    currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,dRDir{dRDirIndx}];
                    destDir=[destinationRoot,filesep,rDDir{rDDirIndx},filesep,dRDir{dRDirIndx}];
                    
                    mkdir(destDir);
                    
                    
                    
                    
                    %The top level directory is that of receptor density
                    
                    %copy rates and densities file to new directory
                    %append file name with number indicating movie #
                    copyfile(currDir,destDir);
                    
                    
                    
                   
%                 end
            end
        end
    end
end