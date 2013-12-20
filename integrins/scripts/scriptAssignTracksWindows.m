
for i = 1 : length(movieStructAlphaVY773A)
    
    disp(num2str(i))
    
    if movieStructAlphaVY773A(i).activityLevel > 1
        
        topDir = movieStructAlphaVY773A(i).fileName{1};
        %         topDir = topDir(11:end);
        %         cd([topDir '/analysisAlphaVY773A/furtherAnalysis'])
        tmp = regexprep(topDir,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        cd([topDir '\analysisAlphaVY773A\furtherAnalysis'])
        %         %%% for randomization test
        %         cd([topDir '\analysisAlphaVY773A\furtherAnalysis\randomizationTest'])
        %         %%% for randomization test
        
        mkdir adaptiveWindowsSym
        cd adaptiveWindowsSym
        
        lengthMinMax = [5 99];
        
        %         firstWindowFile = [topDir '/analysisCellEdgeSmall/windows/windows_frame__frame_01.mat'];
        firstWindowFile = [topDir '\analysisCellEdgeSmall\windows\windows_frame__frame_01.mat'];
        windowsAll = putWindowsTogether(firstWindowFile);
        load ../../../analysisCellEdgeSmall/protrusion_samples/protrusion_samples.mat
        %         %%% for randomization test
        %         load ../../../../analysisCellEdgeSmall/protrusion_samples/protrusion_samples.mat
        %         %%% for randomization test
        
        load ../tracksDiffusionLength5InMask.mat
        %         %%% for randomization test
        %         load ../tracksLength5InMaskRandom.mat
        %         %%% for randomization test
        
        maxFrame = (size(windowsAll,1)-1)*400;
        
        [windowTrackAssign,trackWindowAssign,trackWindowAssignComp,windowTrackAssignExt] = ...
            assignTracks2Windows(tracksFinal,windowsAll,1:400:maxFrame+1,1);
        
        %         if i ~= 16
        
        save('windowsActivityTracks','protSamples','windowsAll','trackWindowAssign',...
            'trackWindowAssignComp','windowTrackAssign','windowTrackAssignExt')
        
        %         else
        %
        %             save('windowsActivityTracks1','protSamples','windowsAll')
        %             save('windowsActivityTracks2','trackWindowAssign','trackWindowAssignComp','windowTrackAssign')
        %             size1 = size(windowTrackAssignExt,4);
        %             size1 = floor(size1/2);
        %             windowTrackAssignExt1 = windowTrackAssignExt(:,:,:,1:size1);
        %             windowTrackAssignExt2 = windowTrackAssignExt(:,:,:,size1+1:end);
        %             save('windowsActivityTracks3','windowTrackAssignExt1')
        %             save('windowsActivityTracks4','windowTrackAssignExt2')
        %             clear windowTrackAssignExt1 windowTrackAssignExt2
        %
        %         end
        
    end
    
end
