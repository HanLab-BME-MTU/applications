function [ output_args ] = GCAmakeInvertedRawWrapper(projList)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%  for iGroup = 1:numel(toPlot.info.names);
projList 
%     projList = toPlot.info.projList{iGroup};
   nProjs = size(projList,1); 
    for iProj = 1:nProjs
        load([projList{iProj,1} filesep 'ANALYSIS' filesep  'movieData.mat']);
        imDir = MD.getChannelPaths{1};
        saveDir = [MD.outputDirectory_  filesep  'Visualization_Overlays' filesep 'Raw'];
        
        if ~isdir(saveDir)
            mkdir(saveDir)
            run = 1; 
        else 
            run =0; 
        end
        if run == 1;
      makeNiceMovie(imDir,saveDir,0.216,16,0,[]);
      %% use ffmpeg to make movie
        %cd(['C:\Users\Maria\Desktop\' ...
         %   'ffmpeg-20141210-git-ae81680-win64-static\ffmpeg-20141210-git-ae81680-win64-static\bin']);
        
        %name = projList{iProj,2};
      %  ffmpeg -i pathToFrames/frame_%0Nd.tif -r 15 -b 20000k filename.mp4
        %execute = ['ffmpeg -r 5 -i ' saveDir filesep '%03d.png' ...
        %' -crf 22 -pix_fmt yuv420p -b 20000k ' saveDir filesep name '.mp4'];
        %system(execute)
        %cd(saveDir)
        else 
            display(['Image Directory Already Found for ' MD.outputDirectory_ ' Skipping']); 
    end
    
% end
end

