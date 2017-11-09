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
            R = ColorStackReader(path);
            for i=1:R.getSizeC;
                C(i) = Channel(path);
            end;
            if(isdir(path))
                outputDirectory = path;
            else
                outputDirectory = fileparts(path);
            end
            obj = obj@MovieData(C,outputDirectory);
        end
        function R = initReader(obj)
            R = ColorStackReader(obj.channels_(1).channelPath_);
        end
        function b = isBF(obj)
            b = true;
        end
    end
    methods ( Static )
        function ML = getMovieList(path)
            D = dir([path filesep '*.tif']);
            D2 = D(~[D.isdir]);
            for i = 1:length(D2)
                MD(i) = ColorStackMovieData([path filesep D2(i).name]);
                MD(i).setPath(pwd);
                MD(i).setFilename(['colorStackMovieData' num2str(i) '.mat']);
                MD(i).save;
            end
            ML = MovieList(MD,path);
        end
    end
    
end

