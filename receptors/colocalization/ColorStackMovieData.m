classdef ColorStackMovieData < MovieData
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = ColorStackMovieData(path)
            if(nargin < 1)
                path ='/project/biophysics/jaqaman_lab/vegf_tsp1/touretLab/140723_ActinCD36Fyn/+PAO+TSP/imagesOriginal';
            end
            for i=1:3;
                C(i) = Channel(path);
            end;
            obj = obj@MovieData(C,path);
        end
        function R = initReader(obj)
            R = ColorStackReader(obj.channels_(1).channelPath_);
        end
        function b = isBF(obj)
            b = true;
        end
    end
    
end

