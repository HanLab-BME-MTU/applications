classdef LaminsImageMultiLength < lamins.classes.LaminsImage
    %LAMINSIMAGEMULTILENGTH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        steerableGabor
        orientationExtendedSkel
    end
    
    methods
        function obj = LaminsImageMultiLength(varargin)
            if(isa(varargin{1},'lamins.classes.LaminsImage'))
                L = varargin{1};
                varargin = { L.parent, L.coordinates{:} };
            end
            obj = obj@lamins.classes.LaminsImage(varargin{:});
        end
        function g = get.steerableGabor(obj)
            sigma = 5;
            nAngles = 36;
            if(isempty(obj.steerableGabor))
                g(3) = struct('response',[],'nms',[],'theta',[],'a',[],'l',[]);
                for i=1:3
                    [g(i).response,g(i).nms,g(i).theta,g(i).a,g(i).l] = steerableGaborFilter(double(obj),sigma,i,[],nAngles);
                end
                obj.steerableGabor = g;
            else
                g = obj.steerableGabor;
            end
        end
        function eskel = get.orientationExtendedSkel(obj)
            import lamins.functions.*;
            if(isempty(obj.orientationExtendedSkel))
                endpts = bwmorph(obj.threshSkel,'endpoints');
                vector = getEndPointOrientationVectors(obj.threshSkel,obj.steerableGabor);
                maxPixelExtension = 20;
                connectivity = 8;
                obj.extendedSkel = extendVectorUntilConnected(obj.threshSkel,endpts,vector,maxPixelExtension,connectivity);
            end
            eskel = obj.orientationExtendedSkel;
        end
    end
    
end

