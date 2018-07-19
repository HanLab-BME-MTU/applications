classdef LaminsData < matlab.mixin.SetGet
    %LaminsData Represents a movieData containing information about Lamins
    %A, B1, and B2
    %   Detailed explanation goes here
    
    
    properties
        movieData;
        cellReader;
        params;
        outputDirectory_;
        thumbs;
        reader;
    end
    
    methods
        function obj= LaminsData(MD)
            import lamins.functions.*;
            if(nargin > 0)
                obj.movieData = MD;
                obj.cellReader = CellReader(CachedReader(MD.getReader()));
                obj.params = laminParams(MD);
                obj.outputDirectory_ =  [MD.outputDirectory_ filesep 'lamins'];
            end
        end
        function R = get.reader(obj)
            R = obj.cellReader;
        end
        corrmatrix = calcCorrMatrix(obj,cache)
        plotCorrMatrix(obj)
        thumbs = getThumbs(obj,scale,subs)
        showThumbs(obj)
        S = getSkeleton(obj,varargin)
        [edges_cc,vertices_cc, pairs, skel] = getSkeletonGraph(obj,varargin)
        showSkeleton(obj,varargin)
        [mask,thresh] = getNuclearMask(obj,varargin)
        plotMaskedCdf(obj,varargin)
        images = getImages(obj)
    end
end
