
%Khuloud Jaqaman, May 2015

sourceRoot = '/project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/09/091814/probeISruns';
saveSource='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170411/probe';
%Define strings for directory hierarchy as needed
rDDir = {'rD4'};%,'rD10','rD16'
aPDir = {'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};
outDirNum = 6:10;
lRDir = {'lR1p0'};

%Define number of label ratio
numLabelRatio = length(lRDir);

fprintf('\n===============================================================');

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
    tic    
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
        
        fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
        
        compTracksVec = cell(numLabelRatio,1);
        
        %Original output is not organized by label ratio since
        %receptorInfoLabeled for each label ratio is saved as a struct.
        %Iterate through the outputs, to pull out each receptorInfoLabeled
        %then the compTracks.  
        for outDirIndx = 1 : length(outDirNum)
            
            %name of current directory
            currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx))];
            
            %Load receptorInfoLabeled
            tempRecepInfo = load([currDir,filesep,...
                'receptorInfoAll',int2str(outDirNum(outDirIndx)),'.mat']);
            
            %generate receptorInfoLabeled for lR = 1
            receptorInfoLabeled = genReceptorInfoLabeled(tempRecepInfo.receptorInfoAll,1,[1 0.3]);
            
            %save
            currDir1 = [saveSource,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx))];
            save([currDir1,['/receptorInfoLabeledExtra' int2str(outDirNum(outDirIndx))]],'receptorInfoLabeled','-v7.3');
             mkdir(curDir1)
             
            %Pull out compTracks
            [compTracksVec{:}] = receptorInfoLabeled(1:numLabelRatio).compTracks;

            %For each label ratio, the inner most directory, create the
            %directory and save compTracks.            
            for lRDirIndx=1:numLabelRatio
                
                fprintf('\n   Out = %d, lR = %s ',outDirIndx,lRDir{lRDirIndx});
                
                currOutDir = [currDir1,filesep,lRDir{lRDirIndx}];
                
                %Create the direcotry
                mkdir(currOutDir)
                
                %Write compTracks
                compTracks = compTracksVec{lRDirIndx};
                save([currOutDir,'/compTracks'],'compTracks','-v7.3');
                
                fprintf('... done.');
                
            end %for each labelRatio
            
            clear compTracks tempRecepInfo receptorInfoLabeled
            
        end %for each outDir
        
        clear compTracksVec
                            
    end %for each aP
    
    elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    
end %for each rD

fprintf('\n\nAll done.\n');

clear



    
