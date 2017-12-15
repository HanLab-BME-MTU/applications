classdef Animation < hgsetget & matlab.mixin.Copyable & handle
%% This class encapsulate a 2D dynamic or static animation, it can be:
%% - A MIP
%% - An overlay
%% - Graph
    methods (Abstract)
        img=loadView(obj,fIdx)
        n=getFrameNb(obj)

    end

    
    methods
        function saveVideo(obj,pathToVideoFile)
            video = VideoWriter(pathToVideoFile);
            video.FrameRate = 5;  % Default 30
            video.Quality = 100;    % Default 75
            open(video)
            for fIdx=1:obj.getFrameNb();
                img=obj.loadView(fIdx);
                writeVideo(video,img);
            end
            close(video)
        end


    end
end