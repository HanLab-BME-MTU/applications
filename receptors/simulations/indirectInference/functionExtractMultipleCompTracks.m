function functionExtractMultipleCompTracks(sourceRoot,rDDir,aPDir,dRDir,outDirNum,lRDir)

%function to extract compTracks from receptorInfoLabeled located in the
%indicated source folders. Note that receptorInfoLabeled has multiple entries
%for different label ratios used, which will each get its own directory.
%
%SYNOPSIS functionExtractMultipleCompTracksdR(sourceRoot,rDDir,aPDir,dRDir,outDirNum,lRDir)
%
%INPUT
%  sourceRoot : path for the receptorInfoLabeled.
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
% OUTPUT
% There is no output, the results will be saved in the sourceRoot directory
% as compTracks.m
%
% Luciana de Oliveira, April 2017


%Define number of label ratio
numLabelRatio = length(lRDir);

fprintf('\n===============================================================');

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
    tic
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
        
        
        for dRDirIndx = 1 : length(dRDir)
            
            
%             fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
            
            compTracksVec = cell(numLabelRatio,1);
            
            %Original output is not organized by label ratio since
            %receptorInfoLabeled for each label ratio is saved as a struct.
            %Iterate through the outputs, to pull out each receptorInfoLabeled
            %then the compTracks.
            for outDirIndx = 1 : length(outDirNum)
                
                %name of current directory
                currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
                    dRDir{dRDirIndx},filesep, aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx))];
                
                %Load receptorInfoLabeled
                
                tempRecepInfo = load([currDir,filesep,...
                    'receptorInfoLabeled',int2str(outDirNum(outDirIndx)),'.mat']);
                
                %Pull out compTracks for each labelRatio defined above
                [compTracksVec{:}] = tempRecepInfo.receptorInfoLabeled(1:numLabelRatio).compTracks;
                
                %For each label ratio, the inner most directory, create the
                %directory and save compTracks.
                for lRDirIndx=1:numLabelRatio
                    
                    fprintf('\n   Out = %d, lR = %s ',outDirIndx,lRDir{lRDirIndx});
                    
                    currOutDir = [currDir,filesep,lRDir{lRDirIndx}];
                    
                    %Create the direcotry
%                     mkdir(currOutDir)
                    
                    %Write compTracks
                    compTracks = compTracksVec{lRDirIndx};
                    save([currOutDir,'/compTracks'],'compTracks','-v7.3');
                    
                    clear compTracks
                    
                    fprintf('... done.');
                    
                end %for each labelRatio
            end
            
            clear compTracks tempRecepInfo
            
        end %for each outDir
        
        clear compTracksVec
        
    end %for each aP
    
    elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    
end %for each rD

fprintf('\n\nAll done.\n');

clear




