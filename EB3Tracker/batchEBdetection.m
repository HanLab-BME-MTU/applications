function batchEBdetection(doDetect,doTrack)

close all

[fileName,pathName] = uigetfile('*.mat','Select projList.mat file containing projects to run');
load([pathName filesep fileName]);
clear 'fileName' 'pathName'

for i=1:size(projList,1)
    try

        
        if doDetect==1
            [movieInfo]=eb1SpotDetector(projList(i),[],16,1,1);
        end
        if doTrack==1
            homeDir=pwd;
            [tempDir] = formatPath(projList(i).anDir);
            cd([tempDir filesep 'feat'])
            load 'movieInfo.mat'
            cd ..
            scriptTrack_EB3 
            cd(homeDir)
        end
        
         
        

       



    catch
        disp('problem with')
        projList(i).anDir
    
    end
end

close all







