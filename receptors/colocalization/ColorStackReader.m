classdef ColorStackReader < Reader
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        path;
    end
    
    methods
        function obj = ColorStackReader(path, varargin)
            obj.path = path;
            if isempty(varargin) 
                fileNames = obj.getImageFileNames();
            else
                fileNames = obj.getImageFileNames(varargin{:});
            end
            
            imageInfo = imfinfo([path filesep fileNames{1}]);
            obj.sizeC = length(imageInfo);
            obj.sizeX = imageInfo(1).Width;
            obj.sizeY = imageInfo(1).Height;
            obj.sizeZ = length(fileNames);
            obj.bitDepth = imageInfo(1).BitDepth;
            obj.sizeT = 1;
           
        end
        function sizeX = getSizeX(obj,varargin)
            sizeX = obj.sizeX;
        end
        function sizeY = getSizeY(obj,varargin)
            sizeY = obj.sizeY;
        end
        function sizeZ = getSizeZ(obj,varargin) 
            sizeZ = obj.sizeZ;
        end
        function sizeC = getSizeC(obj,varargin)
            sizeC = obj.sizeC;
        end
        function sizeT = getSizeT(obj,varargin)
            sizeT = obj.sizeT;
        end
        function bitDepth = getBitDepth(obj,varargin)
            bitDepth = obj.bitDepth;
        end
        function imageFileNames = getImageFileNames(obj,z) %Changes
            D = dir([obj.path filesep '*.tif']);
            D2 = D(~[D.isdir]);
            if(nargin < 2)
                z = 1:length(D2);
            end
            imageFileNames = { D2(z).name };
        end
        function name = getChannelNames(obj,c)
            names = { '1' '2' '3'};
            name = names(c);
%             name = {num2str(c)};
        end
        function I = loadImage(obj,c,t,z) %Changes
            assert(t == 1,'t is out of range');
            if(nargin < 3)
                t = 1;
            end
            if(nargin < 4)
                z = 1;
            end
            filename = obj.getImageFileNames(z);
            I = imread([obj.path filesep filename{1}],c);
        end
    end
    
end

