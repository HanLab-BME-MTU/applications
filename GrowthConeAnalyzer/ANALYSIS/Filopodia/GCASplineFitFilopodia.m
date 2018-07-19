function [ filoInfo ] = GCASplineFitFilopodia(filoInfo,varargin)
%GCASplineFitFilopodia

%%
%%  INPUT:
%    filoInfo: (REQUIRED) : rx1 structure
%    where r is the number of 'filopodia' detected per frame
%    Note: Each linear piece of a thin branch is likewise considered a 'filopodia' in this
%    in this structure. It is distinguished from veil by a field designation
%    'filopodia' are grouped into branches by another field ID.
%
%
%
% PARAMS:
% 'SplineTolerance' (PARAM) :  Positive scalar
%     The tolerance value to use when fitting a spline to filopodia coordinates
%     Larger numbers will smooth the filo more, while 0 will use an interpolating spline.
%     See spaps for details.
%     Default is 2
%     (Note the default was chosen when optimizing the
%     curvature calcs empirically, larger values than this tended to smooth
%     out the subtle curvature that was apparent visually while smaller
%     values picked up on some of the curvature due to pixelation).
%     (Note MB: check one more time before release)
%
%
% OUTPUT:
% filoInfo:  : rx1 structure with new field
%  with
%  ... OLD fields
%  ...
%  and NEW Field
%  .Ext_coordsXY_SplineFit  = rx2 double array
%     where r (rows) is the number of coordinates in the filo
%     and c (cols) is the vector of corresponding coords (x,y)
%     corresponding to the spline fit results
%
%% Check Input
ip= inputParser;

ip.CaseSensitive = false;
%REQUIRED

ip.addRequired('filoInfo');

%
ip.addParameter('SplineTolerance',2); %
ip.addParameter('TSOverlays',false); 
ip.parse(filoInfo,varargin{:});

%% START

for iFilo = 1:length(filoInfo)
    
    % truncate at the fit point
    pixIndices = filoInfo(iFilo).('Ext_pixIndices');
    
    idxEnd = find(pixIndices == filoInfo(iFilo).([ 'Ext_endpointCoordFitPix']));
    vertices = filoInfo(iFilo).('Ext_coordsXY');
    vertices = vertices(1:idxEnd,:);
    
    if ip.Results.TSOverlays 
       scatter(vertices(:,1),vertices(:,2),1,'r','filled'); 
    end 
    
    if size(vertices,1)>2 % will not fit is the length of the points is < 2
        
        [~,vertices]   = spaps( 1:size(vertices,1), vertices',ip.Results.SplineTolerance);
        vertices = vertices';
        filoInfo(iFilo).Ext_coordsXY_SplineFit = vertices;
        if ip.Results.TSOverlays 
            scatter(vertices(:,1),vertices(:,2),1,'k'); 
            %plot(vertices(:,1),vertices(:,2),'color','k'); 
        end 
    else
        filoInfo(iFilo).Ext_coordsXY_SplineFit = [NaN NaN];
    end
    
end



end







