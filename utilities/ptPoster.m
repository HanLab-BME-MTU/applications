function ptPoster (hObject)
% ptPoster       gathers some information for the poster we are making
%
% SYNOPSIS       ptPoster (hObject)
%
% INPUT          none
%                
% OUTPUT         useful mat files          
%
% DEPENDENCIES   ptPoster  uses {nothing}
%                                  
%                ptPoster is used by { PolyTrack_PP }
%
% Andre Kerstens, Mar 04
handles=guidata(hObject);

% Starting frame relativ to first frame analysed with polytrack
startFrame = round((handles.postpro.analfirstimg - handles.jobvalues.firstimage) / handles.jobvalues.increment) + 1;

% Make sure we don't get negative frame numbers
if startFrame < 1
    startFrame = 1;
end

% Last frame relative to first frame analysed with polytrack
stopFrame = floor((handles.postpro.anallastimg - handles.jobvalues.firstimage) / handles.jobvalues.increment + 0.00001) + 1;

% Make sure we stop at the last image and not further than that one
if stopFrame > handles.jobvalues.lastimage
    stopFrame = handles.jobvalues.lastimage;
end

numberOfFrames = stopFrame - startFrame + 1;

saveAllPath=handles.postpro.saveallpath;

% Load all needed values
load ('percentageSingle');
load ('averageDisplacement');
load ('displacementVel');
load ('averageDisplacementVel');
load ('averageSingleDisplacement');
load ('averageClusterDisplacement');
load ('displacement');
load ('velocityVariance');
%load ('thresholdedCells');

% Generate a plot with both average velocity and variance of single cells and 
% cells within clusters in one figure
% h_fig = figure;
% subplot (2,1,1); plot (averageClusterDisplacement), title ('Average velocity of cells within clusters');
% xlabel ('Frames');
% ylabel ('displacement per frame (in pixel)');
% ymax = max (averageClusterDisplacement);
% axis ([startFrame stopFrame 0 ymax]);
% subplot (2,1,2); plot (averageSingleDisplacement), title('Average velocity of single cells');
% xlabel ('Frames');
% ylabel ('displacement per frame (in pixel)');
% ymax = max(averageSingleDisplacement);
% axis ([startFrame stopFrame 0 ymax]);

% hgsave (h_fig, [saveAllPath filesep 'velocitySingleCluster.fig']);
% print (h_fig, [saveAllPath filesep 'velocitySingleCluster.eps'], '-depsc2', '-tiff');
% print (h_fig, [saveAllPath filesep 'velocitySingleCluster.tif'], '-dtiff');

% average the average velocity of cells
winWidth = 25;
winSigma = 0.5;
window = gausswin (winWidth, winSigma);
averageDisplacementVelFiltered = filtfilt (window ./ sum (window), 1, averageDisplacementVel);
averageDisplacementFiltered = filtfilt (window ./ sum (window), 1, averageDisplacement);
velocityVarianceFiltered = filtfilt (window ./ sum (window), 1, velocityVariance);

% Calculate start, end and cutover values for the average velocity (start
% and end are averaged over 20 values)
avgDispVelStart = sum (averageDisplacementVel(1:20)) / 20;
avgDispVelEnd = sum (averageDisplacementVel(end-20:end)) / 20;
avgDispVelCutover = (abs (avgDispVelEnd - avgDispVelStart) / 2) + min ([avgDispVelStart avgDispVelEnd]);
avgDispVelCutoverFrame = find (abs(averageDisplacementVelFiltered - avgDispVelCutover) == min (abs (averageDisplacementVelFiltered - avgDispVelCutover)));

% Calculate start, end and cutover values for the three frame displacement (start
% and end are averaged over 20 values)
avgDispStart = sum (averageDisplacement(1:20)) / 20;
avgDispEnd = sum (averageDisplacement(end-20:end)) / 20;
avgDispCutover = (abs (avgDispEnd - avgDispStart) / 2) + min ([avgDispStart avgDispEnd]);
avgDispCutoverFrame = find (abs(averageDisplacementFiltered - avgDispCutover) == min (abs (averageDisplacementFiltered - avgDispCutover)));

avgDispValues = [avgDispVelStart avgDispVelEnd avgDispVelCutover avgDispVelCutoverFrame; avgDispStart avgDispEnd avgDispCutover avgDispCutoverFrame];

% Prepare vector that holds single and cluster cell numbers measured
% against thresholds
thresholdedCells = zeros(numberOfFrames-1,4);

for iFrame = startFrame : (stopFrame - 1)    
   % Now we'll calculate single and cluster cell numbers via another method
   realDisplacement = displacement (find (displacement (:,iFrame)), iFrame);
   %thresholdValue= averageDisplacementFiltered (iFrame,1);
   thresholdValue = avgDispEnd;
   %thresholdValue = 10;
   thresholdedCells (iFrame,1) = sum (realDisplacement > thresholdValue);  % Store nr of single cells
   thresholdedCells (iFrame,2) = sum (realDisplacement <= thresholdValue); % Store nr of cluster cells
   thresholdedCells (iFrame,3) = thresholdedCells (iFrame,1) / length (realDisplacement); % perc. single cells
   thresholdedCells (iFrame,4) = thresholdedCells (iFrame,2) / length (realDisplacement); % perc. cluster cells
end

% Let's average this stuff again
winWidth = 25;
winSigma = 0.5;
window = gausswin (winWidth, winSigma);
thresholdedCellsFiltered = filtfilt (window ./ sum (window), 1, thresholdedCells);

% Generate a plot with both average velocity and variance of cells in one figure
h_fig = figure;
%subplot (2,1,1); plot (averageDisplacement), title('Average velocity of cells');
plot (averageDisplacementVel), title('Average velocity of cells');
hold on;
plot (averageDisplacementVelFiltered, 'r');
plot (avgDispVelCutover, 'k');
xlabel ('Frames');
ylabel ('displacement per frame (in pixel)');
ymax = max(averageDisplacementVel);
axis ([startFrame stopFrame 0 ymax]);
hold off;
% subplot (2,1,2); plot (velocityVariance), title('Variance of avarage velocity of cells');
% hold on;
% plot (velocityVarianceFiltered, 'r');
% xlabel ('Frames');
% ylabel ('variance');
% ymax = max(velocityVariance);
% axis ([startFrame stopFrame 0 ymax]);
% hold off;

% hgsave (h_fig, [saveAllPath filesep 'velocityVariance.fig']);
% print (h_fig, [saveAllPath filesep 'velocityVariance.eps'], '-depsc2', '-tiff');
% print (h_fig, [saveAllPath filesep 'velocityVariance.tif'], '-dtiff');

% average the percentage of single cells
winWidth = 25;
winSigma = 0.5;
window = gausswin (winWidth, winSigma);
percentageSingleFiltered = filtfilt (window ./ sum (window), 1, percentageSingle);

% Calculate start, end and cutover values for percentage single cells (start and end are averaged over
% 20 values)
percSingleStart = sum (percentageSingle(1:20)) / 20;
percSingleEnd = sum (percentageSingle(end-20:end)) / 20;
percSingleCutover = (abs (percSingleEnd - percSingleStart) / 2) + min ([percSingleStart percSingleEnd]);
percSingleCutoverFrame = find (abs(percentageSingleFiltered - percSingleCutover) == min (abs (percentageSingleFiltered - percSingleCutover)));
percSingleValues = [percSingleStart percSingleEnd percSingleCutover percSingleCutoverFrame];

% Plot the percentage of single cells calculated by colin's method
h_fig=figure,title(handles.jobvalues.imagename);
ymax=max(percentageSingle);
subplot(1,1,1); plot(percentageSingle), title('Percentage of single cells (colins method)');
hold on;
plot(percentageSingleFiltered, 'r');
xlabel('Frames');
ylabel('percentage of single cells');
axis([0 stopFrame 0 1]);
hold off;
           
hgsave(h_fig,[saveAllPath filesep 'percentSingle.fig']);
print(h_fig, [saveAllPath filesep 'percentSingle.eps'],'-depsc2','-tiff');
print(h_fig, [saveAllPath filesep 'percentSingle.tif'],'-dtiff');
	
% Plot all the cells in all the frames against their velocities and put the
% filtered average velocity there as well
transposedDisplacementVel = displacementVel';
h_fig = figure;
plot (transposedDisplacementVel,'rx'), title('Velocity of cells (black line shows average velocity)');
hold on;
plot (averageDisplacementVelFiltered, 'k');
xlabel ('Frames');
ylabel ('displacement per frame (in pixel)');
hold off;

% Plot all the cells in all the frames against their three frame displacements and put the
% filtered average displacement there as well
transposedDisplacement = displacement';
h_fig = figure;
plot (transposedDisplacement,'rx'), title('three frame displacements (black line shows average)');
hold on;
plot (averageDisplacementFiltered, 'k');
xlabel ('Frames');
ylabel ('displacement per three frames (in pixel)');
hold off;

% Plot number of single cells with new method
h_fig = figure;
% ymax=max(thresholdedCells(:,1));
% subplot (2,1,1); plot (thresholdedCells(:,1),'g'), title('Number of single cells');
% hold on;
% plot (thresholdedCellsFiltered(:,1),'r');
% xlabel ('Frames');
% ylabel ('# cells');
% hold off;
% Plot perc. of single cells (calculated with vel. method)
%subplot (2,1,2); plot (thresholdedCells(:,3),'g'), title('Percentage of single cells');
plot (thresholdedCells(:,3),'g'), title('Percentage of single cells');
hold on;
plot (thresholdedCellsFiltered(:,3),'r');
ymax=max(thresholdedCells(:,3));
xlabel ('Frames');
ylabel ('Percentage');
hold off;

% Plot number of cluster cells
% h_fig = figure;
% ymax=max(thresholdedCells(:,2));
% subplot (2,1,1); plot (thresholdedCells(:,2),'g'), title('Number of cluster cells');
% hold on;
% plot (thresholdedCellsFiltered(:,2),'r');
% xlabel ('Frames');
% ylabel ('# cells');
% hold off;
% % Plot perc. of cluster cells (calculated with vel. method)
% subplot (2,1,2); plot (thresholdedCells(:,4),'g'), title('Percentage of cluster cells');
% hold on;
% plot (thresholdedCellsFiltered(:,4),'r');
% ymax=max(thresholdedCells(:,4));
% xlabel ('Frames');
% ylabel ('Percentage');
% hold off;

guidata(hObject, handles);

% Write values to the disk as mat file
save ('transposedDisplacementVel','transposedDisplacementVel');
save ('averageDisplacementVelFiltered','averageDisplacementVelFiltered');
save ('averageDisplacementFiltered','averageDisplacementFiltered');
save ('velocityVarianceFiltered','velocityVarianceFiltered');
save ('percentageSingleFiltered','percentageSingleFiltered');
save ('percSingleValues','percSingleValues');
save ('avgDispValues','avgDispValues');

% Write the thresholdedCells matrix, the avg velocity and percent single cell values 
% to the disk in ascii format
fid = fopen ('thresholdedCells.txt','w');
fprintf (fid, '%f %f %f %f\n', thresholdedCells');
fclose (fid);

fid = fopen ('percSingleValues.txt','w');
fprintf (fid, '%f %f %f %f\n', percSingleValues');
fclose (fid);

fid = fopen ('avgDispValues.txt','w');
fprintf (fid, '%f %f %f %f %f %f %f %f\n', avgDispValues');
fclose (fid);
