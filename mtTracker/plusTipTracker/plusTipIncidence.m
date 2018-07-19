function projData=plusTipIncidence(projData,xCoord,yCoord)
% plusTipIncidence read the project tracks and the cell boundary, get the incident angles of the mt with the cell boundary and plot them
%
% Input:
%   projData  : structure containing dir and cell boundary data
%   xCoord, yCoord: the final tracks in current ROI
% Output:
%   projData: Data added with incident angle results

% Liya Ding, Jan 2012

% get the size and the tracks
[nTracks nFrames] = size(xCoord(:,1));

% Initialize the vector and arrays to be used
normalVectors = zeros(nTracks,2);
mtDirectionVectors = zeros(nTracks,2);
incidence = zeros(nTracks,1);

% Get the contour and normal calculated from projData
contourYX = projData.contourYX;
normalYX = projData.normalYX;

% Get the roiMask but redraw with the cell boundary info.
poi_x_dir = getFileNameBody(getFileNameBody(projData.anDir));
roiMask = imread([poi_x_dir filesep 'roiMask.tif']);
roiMask = roipoly(roiMask,contourYX(:,2),contourYX(:,1));
dist_map = bwdist(roiMask) + bwdist(1-roiMask)-1;

display_flag=1;

if display_flag ==1
h = figure; axis image;
end

% For every track
for iTrack = 1:nTracks
    thistrack_x = xCoord(iTrack,:);
    thistrack_y = yCoord(iTrack,:);
    
    thistrack_x = round(thistrack_x(~isnan(thistrack_y)));
    thistrack_y = round(thistrack_y(~isnan(thistrack_y)));
    
    if display_flag ==1
    hold on;
    plot(thistrack_x, thistrack_y,'m');
    end
    
    ind = sub2ind(size(dist_map),thistrack_y,thistrack_x);
    
    dist_array = dist_map(ind);
    
    ind_close = zeros(1,2);
    ind_close_1 = find(dist_array == min(dist_array));
    ind_close(1) = ind_close_1(1);
    dist_array(ind_close(1))=Inf;
    
    ind_close_2 = find(dist_array == min(dist_array));
    ind_close(2) = ind_close_2(1);
    
    if(ind_close(2)<ind_close(1))
        ind_close = fliplr(ind_close);
    end
    
    if(ind_close(2)~=1)
    ind_close(1) = ind_close(2)-1;
    end
    
    % Get the two points closest to the cell boundary
    close_x = [thistrack_x(ind_close(1)) thistrack_x(ind_close(2))];
    close_y = [thistrack_y(ind_close(1)) thistrack_y(ind_close(2))];
    
    % the direction vector of the mt movement
    mtDirection = [close_y(2)-close_y(1) close_x(2)-close_x(1)];
    mtDirectionVectors(iTrack,1:2) = mtDirection/norm(mtDirection);
    for_display_mt(iTrack,1:2) = [close_y(2) close_x(2)];
        
    % normal of the cell boundary on the point closest to mt points
    ind_normal = findthepoint(close_x,close_y,contourYX);
    normalVectors(iTrack,1:2) = normalYX(ind_normal,:);
    for_display_normal(iTrack,1:2) = [contourYX(ind_normal,2) contourYX(ind_normal,1)];
    
    % get the incidence angle from the dot product of the two normalized direction
    % vector
    incidence(iTrack) = acos(dot(mtDirectionVectors(iTrack,1:2),normalVectors(iTrack,1:2)));
end

if display_flag ==1

quiver(for_display_mt(:,2),for_display_mt(:,1),mtDirectionVectors(:,2),mtDirectionVectors(:,1),'b');
hold on;
quiver(for_display_normal(:,1),for_display_normal(:,2),normalVectors(:,2),normalVectors(:,1),'g');
plot(contourYX(:,2),contourYX(:,1),'r');
axis image; axis ij;
saveas(h,[projData.anDir filesep 'incidence.tif']);
saveas(h,[projData.anDir filesep 'incidence.fig']);
close;

h = figure();
hist(incidence);
saveas(h,[projData.anDir filesep 'hist.tif']);
saveas(h,[projData.anDir filesep 'hist.fig']);
close;


end

end


function ind = findthepoint(x,y,contourYX)
Y = contourYX(:,1) - y(2);
X = contourYX(:,2) - x(2);
Z = X.*X + Y.*Y;
ind = find(Z ==min(Z));
ind = ind(1);
end