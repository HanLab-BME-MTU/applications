function [ success,saveDir] = GCAVisualsInvertedRawMovie(MD)
%GCAVisualsInvertedRawMovie 
% small function to 

% for now just keep success to 1 
success = 1; 
 

imDir = MD.getChannelPaths{1};
saveDir = [MD.outputDirectory_  filesep  'Visualization_Overlays' filesep 'Raw'];

% check if the directory exists : if yes skip 
if ~isdir(saveDir)
    mkdir(saveDir)
    run = 1;
else
    run =0;
    
end

% check to make sure you have the files you need 

if run == 1;
    makeNiceMovie(imDir,saveDir,0.216,16,0,[]);
    %% use ffmpeg to make movie (on windows having some trouble with odd sizes) 
%     cd(['C:\Users\Maria\Desktop\' ...
%        'ffmpeg-20141210-git-ae81680-win64-static\ffmpeg-20141210-git-ae81680-win64-static\bin']);
%     
%     name = projList{iProj,2};
%       ffmpeg -i pathToFrames/frame_%0Nd.tif -r 15 -b 20000k filename.mp4
%     execute = ['ffmpeg -r 5 -i ' saveDir filesep '%03d.png' ...
%     ' -crf 22 -pix_fmt yuv420p -b 20000k ' saveDir filesep name '.mp4'];
%     system(execute)
%     cd(saveDir)
else
    display(['Image Directory Already Found for ' MD.outputDirectory_ ' Skipping']);
end




