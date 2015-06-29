%Script to collect compTracks from receptorInfoLabeled located in the
%indicated source folders. New analysis directory hierarchy will be created
%to hold the compTracks. Note that receptorInfoLabeled has multiple entries
%for different label ratios used.
%
%Robel Yirdaw, November 2014
%


sourceRoot = '/project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/11/112114/targetISruns';
%Use the following to redirect output to a different location than the
%location of this script
analysisRoot = '';
%Define strings for directory hierarchy as needed
rDDir = {'rD4'};
aPDir = {'aP0p5'};
lRDir = {'lR0p2';'lR0p4'};

%numSims is the number of simulation runs or repeats - corresponds to the
%number of out* folders.
numSims = 10;
%Define number of label ratio
numLabelRatio = length(lRDir);

fprintf('\n===============================================================');
%The top level directory is that of receptor density
for rDDirIndx=1:length(rDDir)
    
    tic    
    %Iterate through association probability values per density
    for aPDirIndx=1:length(aPDir)
        fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
        
        compTracksVec = cell(numLabelRatio,1);
        
        %Original output is not organized by label ratio since
        %receptorInfoLabeled for each label ratio is saved as a struct.
        %Iterate through the outputs, to pull out each receptorInfoLabeled
        %then the compTracks.  
        for outDirIndx=1:numSims
            %Load receptorInfoLabeled
            tempRecepInfo = load([sourceRoot,'/',rDDir{rDDirIndx},'/',...
                aPDir{aPDirIndx},'/','out',int2str(outDirIndx),'/',...
                'receptorInfoLabeled',int2str(outDirIndx),'.mat']);
            
            %Pull out compTracks for each labelRatio defined above
            [compTracksVec{:}] = tempRecepInfo.receptorInfoLabeled(1:numLabelRatio).compTracks;

            %For each label ratio, the inner most directory, create the
            %directory and save compTracks.            
            for lRDirIndx=1:numLabelRatio   
                fprintf('\n   Out = %d, lR = %s ',outDirIndx,lRDir{lRDirIndx});
                
                currOutDir = [rDDir{rDDirIndx},'/',...
                    aPDir{aPDirIndx},'/',lRDir{lRDirIndx},'/out',int2str(outDirIndx)];
                
                %Create the direcotry
                mkdir(currOutDir)
                
                %Write compTracks
                compTracks = compTracksVec{lRDirIndx};
                save([currOutDir,'/compTracks',int2str(outDirIndx)],'compTracks','-v7.3');
                
                clear compTracks
                
                fprintf('... done.');
            end %for each labelRatio
            
            clear compTracks tempRecepInfo
            
        end %for each outDir
        
        clear compTracksVec
                            
    end %for each aP
    
    elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    
end %for each rD

fprintf('\n\nAll done.\n');

clear



    
