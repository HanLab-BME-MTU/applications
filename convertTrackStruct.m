function tracks = convertTrackStruct(tracks,skelIn,formatStr)
%CONVERTTRACKSTRUCT converts the input temporal track graph to more intuitive XYZT format
%
% tracksOut = convertTrackStruct(tracksIn,skelIn,formatStr)
%
%       formatStr Options:
%
%           'Imaris' - This gives XYZ coordinates which match with the
%           input temporal edges for plotting in imaris.
%
% Hunter Elliott
% 7-2012
%

%Remove fields for readability
tEdges = tracks.tEdges;
tVerts = tracks.tVerts;

switch formatStr
    
    
    case 'Imaris'
               
        nPts = size(tracks.tVerts,1);
        tracks.XYZ = nan(nPts,3);
        tracks.T = nan(nPts,1);
        
%         nEdge = size(tracks.tEdges,1);
%         
%         for j = 1:nEdge
%            %Starting vertex for this track segment
%            tracks.XYZ(tEdges(j,1),:) = ...
%                skelIn(tVerts(tEdges(j,1),1)).vertices(tVerts(tEdges(j,1),2),:);
%            %Ending vertex for this track segment
%            tracks.XYZ(tEdges(j,2),:) = ...
%                skelIn(tVerts(tEdges(j,2),1)).vertices(tVerts(tEdges(j,2),2),:);                       
%            tracks.T(tEdges(j,1)) = tVerts(tEdges(j,1),1);%Time for starting vertex
%            tracks.T(tEdges(j,2)) = tVerts(tEdges(j,2),1);%Time for end vertex
%         end
        
        for j = 1:nPts
            tracks.XYZ(j,:) = skelIn(tVerts(j,1)).vertices(tVerts(j,2),:);
            tracks.T(j) = tVerts(j,1);
        end
        
%         
%         isNan = isnan(tracks.T);
%         tracks.T = tracks.T(~isNan);
%         tracks.XYZ = tracks.XYZ(~isNan,:);
%                 
    otherwise
        error('unrecognized format string!')
end

