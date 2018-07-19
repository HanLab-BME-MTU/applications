function [ filoInfo ] = GCAAddFilopodiaCurvature(filoInfo)
%GCAAddFilopodiaCurvature

%%
%%  INPUT:
%    filoInfo: (REQUIRED) : rx1 structure
%    where r is the number of 'filopodia' detected per frame
%    Note: Each linear piece of a thin branch is likewise considered a 'filopodia' in this
%    in this structure. It is distinguished from veil by a field designation
%    'filopodia' are grouped into branches by another field ID.
%
% OUTPUT:
% filoInfo:  : rx1 structure with new field
%  with
%  ... OLD fields
%  ...
%  and NEW Field
%  .Ext_FiloCurvature = rx1 double array
%     where r (rows) is the number of coordinates in the filo
%     and c (cols) the local curvature k as calculated by 
%     gcaLineCurvature2D.m (an extern function from )
%% Check Input
ip= inputParser;

ip.CaseSensitive = false;
%REQUIRED

ip.addRequired('filoInfo');

ip.parse(filoInfo);

%% START

for iFilo = 1:length(filoInfo)
    vertices =filoInfo(iFilo).Ext_coordsXY_SplineFit; 
    k = gcaLineCurvature2D(vertices); 
    filoInfo(iFilo).Ext_FiloCurvIndVals = k;   
   
end



end







