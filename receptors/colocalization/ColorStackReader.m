classdef ColorStackReader < Reader
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        path;
        dir;
    end
    
    methods
        function obj = ColorStackReader(path, varargin)
            obj.path = path;
            if(isdir(obj.path))
                obj.dir = obj.path;
            else
                obj.dir = fileparts(obj.path);
            end
            if isempty(varargin) 
                fileNames = obj.getImageFileNames();
            else
                fileNames = obj.getImageFileNames(varargin{:});
            end
            
            imageInfo = imfinfo([obj.dir filesep fileNames{1}{1}]);
            obj.sizeC = length(imageInfo);
            obj.sizeX = imageInfo(1).Width;
            obj.sizeY = imageInfo(1).Height;
            obj.sizeT = length(fileNames{1});
            obj.bitDepth = imageInfo(1).BitDepth;
            obj.sizeZ = 1;
           
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
        function imageFileNames = getImageFileNames(obj,iChan, iFrame, varargin) %Changes
            if(isdir(obj.path))
                D = dir([obj.path filesep '*.tif']);
            else
                D = dir(obj.path);
            end
            D2 = D(~[D.isdir]);
            if(isempty(obj.sizeC))
                if(isdir(obj.path))
                    sizeC = length(imfinfo([obj.path filesep D2(1).name]));
                else
                    sizeC = length(imfinfo(obj.path));
                end
            else
                sizeC = obj.sizeC;          
            end
            imageFileNames = cell(1,sizeC);
            for c = 1:sizeC
                imageFileNames{c} = { D2(:).name };
            end
            if(nargin > 1)
                imageFileNames = imageFileNames{iChan};
            end
            if(nargin > 2)
                imageFileNames = imageFileNames(iFrame);
            end
        end
        function name = getChannelNames(obj,c)
            names = { '1' '2' '3' '4'};
            name = names(c);
%             name = {num2str(c)};
        end
        function I = loadImage(obj,c,t,z) %Changes
            if(nargin < 3)
                t = 1;
            end
            if(nargin < 4)
                z = 1;
            end
            assert(z == 1,'z is out of range');
            filename = obj.getImageFileNames(c,t);
            I = imread([obj.path filesep filename{1}],c);
            if(any(size(I) ~= [obj.sizeY obj.sizeX]))
                I = imresize(I,[obj.sizeY obj.sizeX]);
            end
        end
    end
    
end

