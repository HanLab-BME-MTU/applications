classdef TracksROI < DynROI
    properties (SetAccess = public, GetAccess = public)
    tracks;
    fringe;  
    end

    methods
        function obj = ProjectDynROIProcess(tracks,fringe)
            obj.tracks=tracks;
            obj.fringe=fringe;
        end

        function [minmaxXBorder, minmaxYBorder,minmaxZBorder]=getBoundingBox(obj)
            dynROI=obj.tracks;
            minX=dynROI(1).x;
            minY=dynROI(1).y;
            minZ=dynROI(1).z;
            maxX=0; 
            maxY=0;
            maxZ=0;
            for iP=1:length(dynROI)
                minX=floor(min(min(dynROI(iP).x),minX));
                minY=floor(min(min(dynROI(iP).y),minY));
                minZ=floor(min(min(dynROI(iP).z),minZ));
                maxX=ceil(max(max(dynROI(iP).x),maxX));
                maxY=ceil(max(max(dynROI(iP).y),maxY));
                maxZ=ceil(max(max(dynROI(iP).z),maxZ));
            end
            fringeWidth=obj.fringe;
            maxXBorder=(maxX+fringeWidth);
            maxYBorder=(maxY+fringeWidth);
            maxZBorder=(maxZ+fringeWidth);
            minXBorder=(minX-fringeWidth);
            minYBorder=(minY-fringeWidth);
            minZBorder=(minZ-fringeWidth);

            minmaxXBorder=[minXBorder maxXBorder];
            minmaxYBorder=[minYBorder maxYBorder];
            minmaxZBorder=[minZBorder maxZBorder];
        end

        

    end

end
