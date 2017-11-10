% Initialize variables
emptyFrame = zeros (1,5);

startDir = pwd;
cd (startDir);

posDirs = dir('pos*');
nrPosDirs = length(posDirs);

for movieCount = 1 : nrPosDirs
    
    fprintf (1, '%s - Processing movie: %d\n\n', datestr(now), movieCount);
    
    cd (posDirs(movieCount).name);
    
    resultDirs = dir('result*');
    nrResultDirs = length(resultDirs);
    
    for resultsCount = 1 : nrResultDirs
        
        cd (resultDirs(resultsCount).name);
    
        % Load the cluster (binary) image of the cells (nuclei and halos combined)
        curDir = pwd;
        clusterDir = [curDir filesep 'info'];
        clusterImageDir = [curDir filesep 'body'];
        cd (clusterDir);

        nrFiles = length(dir('clus*'));

        fprintf (1, '%s - Processing frame: ', datestr(now));

        for frameCount = 1 : nrFiles

           fprintf (1, ' %d', frameCount);

           formatStr = sprintf ('%%.%dd', 3);
           imageNr = sprintf (formatStr, frameCount);
           clusterFile = ['clusters' imageNr];
           load (clusterFile);
           
           % Generate figure to plot the convex hull
           %h_fig = figure; imshow (clusterImage);
           %hold on;
           
           % Label the cluster image
           imgLabeledCellArea = bwlabel (clusterImage);

           % Calculate ratio area/convex_hull_area for all cells and clusters
           clear coord;
           regionProps = regionprops (imgLabeledCellArea, 'Area', 'ConvexArea', 'ConvexHull');
           for jCount = 1 : size (regionProps, 1)
              area(jCount) = regionProps(jCount).Area;
              convexArea(jCount) = regionProps(jCount).ConvexArea;
              Solidity(jCount) = area(jCount) / convexArea(jCount);

              %coord = regionProps(jCount).ConvexHull;
              %plot (coord(:,1),coord(:,2),'r-');
           end
           
           % Save and close figure
           %hold off;
           %convexFile = ['clustersWithConvexHull' imageNr '.tif'];
           %print (h_fig, [clusterImageDir filesep convexFile],'-dtiff');
           %close (h_fig);
           
           % Calculate average ratios for the frame
           avgArea = sum(area) / length(area);
           avgConvexArea = sum(convexArea) / length(convexArea);
           avgSolidity = sum(Solidity) / length(Solidity);

           % Properties for the whole frame are stored in frameProps
           % Keep some space for later values as well
           frameProp = zeros (1,5);
           frameProp (1,1) = avgArea;
           frameProp (1,2) = avgConvexArea;
           frameProp (1,3) = avgSolidity;

           % Accumulate the frame properties
           frameProps (1 : size (emptyFrame, 1), 1 : size (emptyFrame, 2), frameCount) = emptyFrame;
           frameProps (1 : size (frameProp, 1), 1 : size (frameProp, 2), frameCount) = frameProp;

        end  % for frameCount

        % Go back to the results dir
        cd (curDir);

        % Save frame properties
        fprintf (1, 'Writing frameProps.mat in %s\n\n', resultDirs(resultsCount).name);
        if exist ('frameProps', 'var')
           save ('frameProps.mat', 'frameProps');
        end
        
        cd ('..');
        
    end  % for resultCount   
    
    cd (startDir);
    
end  % for movieCount