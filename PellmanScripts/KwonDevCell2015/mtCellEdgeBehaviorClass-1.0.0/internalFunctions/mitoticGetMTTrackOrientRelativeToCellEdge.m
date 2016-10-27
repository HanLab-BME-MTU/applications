function [incidence,dispMTVect,dwellMTVect] =mitoticGetMTTrackOrientRelativeToCellEdge(projData,xCoord,yCoord,iMask,varargin) 
% plusTipIncidence read the project tracks and the cell boundary, get the incident angles of the mt with the cell boundary and plot them
%
% Input:
%   projData:          structure containing dir 
%
%   xCoord, yCoord:    the final tracks in current ROI
%   
%   iMask:             the frame number corresponding to the mask you would like to load
%
%   minDispVectOrient: the minimum size of the displacement vector to use
%                      for the orientation calculation (currently in pixels)
%                      the function will search back from the end of the 
%                      MT track until at least this vector length is found
%                      
%   
% Output:
%   [incidence dispMTVect ] : a ntrack x 2 double matrix containing
%   orientation calculations for each tracks and the size ( in pixels ) of
%   the displacement vector used to make that calculation. 

% Liya Ding, Jan 2012
% updated Maria Bagonis August 2013; for pellman collab. 


%% Check Input 
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('minDispVect',3,@(x) isnumeric(x));

ip.addParameter('TSOverlays',true,@(x) islogical(x)); 

ip.parse(varargin{:}); 
    
    
%% INFORMATION 

% get the size and the tracks
[nTracks nFrames] = size(xCoord(:,1));

% Initialize the vector and arrays to be used
normalVectors = zeros(nTracks,2);
mtDirectionVectors = zeros(nTracks,2);
incidence = zeros(nTracks,1);
dispMTVect = zeros(nTracks,1); 
dwellMTVect = zeros(nTracks,1); 

% Get the roiMask but redraw with the cell boundary info.
%poi_x_dir = getFilenameBody(getFilenameBody(projData.anDir));
poi_x_dir  = mitoticUpDirectory(projData.anDir,4); 
if ~isempty(iMask);
    
    masks = mitoticSearchFiles('.tif',[],[poi_x_dir filesep 'masks'],0);
end


if ~isempty(iMask)
    roiMask = logical(imread([masks{iMask,2} filesep  masks{iMask,1}]));
    
else
    
    roiMask = logical(imread([poi_x_dir filesep 'roiMaskSeg.tif']));
end

roiMask = logical(getLargestCC(roiMask));
[contourYX,normalYX] = outercontour_normal(roiMask);
roiMask = roipoly(roiMask,contourYX(:,2),contourYX(:,1));
%dist_map = bwdist(roiMask) + bwdist(1-roiMask)-1;


if ip.Results.TSOverlays
    figure('Visible','off'); axis image;
    imshow(roiMask,[]);
    hold on
end

%% START LOOP 
% For every track
for iTrack = 1:nTracks
    thistrack_x = xCoord(iTrack,:);
    thistrack_y = yCoord(iTrack,:);
    
    %     thistrack_x = round(thistrack_x(~isnan(thistrack_y)));
    %     thistrack_y = round(thistrack_y(~isnan(thistrack_y)));
    % DON't round
    thistrack_x = thistrack_x(~isnan(thistrack_x));
    thistrack_y = thistrack_y(~isnan(thistrack_y));
    
    % plot the track in magenta
    if ip.Results.TSOverlays
        hold on;
        plot(thistrack_x, thistrack_y,'m','LineWidth',1);
        hold on 
        scatter(thistrack_x,thistrack_y,10,'m','filled'); 
    end
   
    ind_close(1) = length(thistrack_x); % let's try just taking the end point

    %% Test for uniform displacement by searching back coordinates
    skipflag=0; % initiate skip flag 
    if ind_close(1) - 1  <= 0 %  if there are not two coordinates you can't take the
        % displacement vector (can happen as if only reading the coordinates
        % int the region ) 
        skipflag =1 ; % set the skip flag to 1 
    else
        %else
        %   back = 1;
        %if ind_close(1) -1 <=0
        back = 1;  % start a single coordinate back
    end
    
    %end
    c = colormap(jet(10));
    if skipflag ~=1
        stopFlag = 0;  % initiate stop flag
        it = 1; % initiate interations
        
        
        while stopFlag == 0
            
            close_x = [thistrack_x(ind_close(1)-back) thistrack_x(ind_close(1))]; % start at first frame back
            close_y = [thistrack_y(ind_close(1)-back) thistrack_y(ind_close(1))];
            
            it = it +1;
            
            % the direction vector of the mt movement
            mtDirection = [close_y(2)-close_y(1) close_x(2)-close_x(1)];
            disp = sqrt((close_y(2) - close_y(1) )^2 + (close_x(2) - close_x(1)) ^2 );
            
            if disp > ip.Results.minDispVect;% 
                if it <=10
                  dispVectPlotColor = c(it,:); 
                else 
                    dispVectPlotColor = c(end,:); % set to highest just to make sure won't error... (want to keep low in general though just to 
                    % make sure have color distinction among points. 
                end 
               
                scatter(close_x,close_y,10,dispVectPlotColor,'filled'); % plot the end points color code them by iteration 
                dispMTVect(iTrack) = disp; % save the displacement used for vector calc
                dwellMTVect(iTrack) = it; % save the number of iterations (the time) 
                stopFlag = 1 ; % disp greater > minDispVect 
            else
                back = back +1 ; % try going a bit farther back so not measuring just noise ( until displacement)
                if ind_close(1) <=0 || ind_close(1)-back <=0
                    skipflag = 1;
                    stopFlag = 1;
                end
            end % end if
        end % end while
    end % if skip
    if skipflag ~=1
        mtDirectionVectors(iTrack,1:2) = mtDirection/norm(mtDirection);
        for_display_mt(iTrack,1:2) = [close_y(2) close_x(2)];
        
        % normal of the cell boundary on the point closest to mt points
        ind_normal = findthepoint(close_x,close_y,contourYX);
        normalVectors(iTrack,1:2) = normalYX(ind_normal,:);
        for_display_normal(iTrack,1:2) = [contourYX(ind_normal,2) contourYX(ind_normal,1)];
        
        % get the incidence angle from the dot product of the two normalized direction
        % vector
        incidence(iTrack) = acos(dot(mtDirectionVectors(iTrack,1:2),normalVectors(iTrack,1:2)));
        
        %
        
        text(thistrack_x(1,1),thistrack_y(1,1),  [(num2str(rad2deg(incidence(iTrack)),2)) ' ' num2str(dispMTVect(iTrack),2)],'color','r','fontSize',16);
        %     if ~isempty(dispValuesFilt)
        %     text(thistrack_x(1,end),thistrack_y(1,end),num2str(dispValuesFilt(iTrack),2),'color','b','fontSize',16);
        %     end
        
    else
        incidence(iTrack)=NaN;
        dispMTVect(iTrack) = NaN;
        dwellMTVect(iTrack) = NaN;
    end % skipFlag
    % plot dispMTVector in Middle of massk
    
end % for iTrack 

%% OPTIONAL PLOTS

if ip.Results.TSOverlays && skipflag ~=1
    
    quiver(for_display_mt(:,2),for_display_mt(:,1),mtDirectionVectors(:,2)*15,mtDirectionVectors(:,1)*15,0,'b','Linewidth',1.5);
    hold on;
    quiver(for_display_normal(:,1),for_display_normal(:,2),normalVectors(:,2)*15,normalVectors(:,1)*15,0,'g','Linewidth',1.5);
    plot(contourYX(:,2),contourYX(:,1),'r');
    axis image; axis ij;
    
    % make a directory 
    cortDirIndiv = [projData.anDir filesep 'meta' filesep 'CorticalInfo' ...
        filesep 'TroubleshootClassifications' filesep 'IncidencePlots' filesep 'minDispVect_' num2str(ip.Results.minDispVect) ];
    if ~isdir(cortDirIndiv)
        mkdir(cortDirIndiv) 
    end 
    
    cortDirCombined = [projData.anDir filesep 'meta' filesep 'CorticalInfo'...
        filesep 'TroubleshootClassifications' filesep 'IncidencePlots' filesep 'All']; 
    if ~isdir(cortDirCombined)
        mkdir(cortDirCombined) 
    end 
    
    fmt = '%03d';
    colormap(gray)
    cd(cortDirIndiv)
    % save mask in individual directory   
    saveas(gcf,[cortDirIndiv filesep 'Mask' num2str(iMask,fmt) '_minDispVectSize' num2str(ip.Results.minDispVect) '.tif']);
    saveas(gcf,[cortDirIndiv filesep 'Mask' num2str(iMask,fmt) '_minDispVectSize' num2str(ip.Results.minDispVect) '.eps'],'psc2'); 
    
    % save mask in combined directory for quick comparison 
    
    saveas(gcf,[cortDirCombined filesep 'Mask' num2str(iMask,fmt) '_minDispVectSize' num2str(ip.Results.minDispVect) '.tif']); 
    close gcf     
end

end


function ind = findthepoint(x,y,contourYX)
Y = contourYX(:,1) - y(2);
X = contourYX(:,2) - x(2);
Z = X.*X + Y.*Y;
ind = find(Z ==min(Z));
ind = ind(1);
end