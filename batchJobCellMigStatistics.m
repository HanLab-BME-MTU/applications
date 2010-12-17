% browse into the folder that contains all the subfolders and execute this
% script.

% set the parameters for the statistics:
display('Now set the parameters for the statistics:');
binPix          =input('Specify the band width for density measurements,  binPix=[  100pix]: ');
minTrackLength  =input('Set the minimal considered track length,  minTrackLength=[12frames]: ');
timeWindow      =input('Time window for averaging the cell velocity, timeWindow =[ 1frames]: ');
dframes         =input('Plot only velocities every ...                          =[24frames]: ');
dt              =input('Specify the time between two frames as fraction of 1h   =[    1/6h]: ');
showTrackerMovie=input('Do you want to plot a tracker movie?        Yes=1 / No=0=[      No]: ');
showEdgeMovie   =input('Do you want to show the tracked sheet edge? Yes=1 / No=0=[      No]: ');

if isempty(binPix)
    binPix=100;
end
if isempty(minTrackLength)
    minTrackLength=12;
end
if isempty(timeWindow)
    timeWindow=1;
end
if isempty(dframes)
    dframes=24;
end
if isempty(dt)
    dt=1/6;
end
if isempty(showTrackerMovie)
    showTrackerMovie=0;
end
if isempty(showEdgeMovie)
    showEdgeMovie=0;
end


% set the folder names:
A = dir();
sizeA = size(A,1);
for k=1:sizeA
    checkVec(k)=(~A(k).isdir || strcmp(A(k).name,'.') || strcmp(A(k).name,'..'));
end
A(checkVec)=[];
% give the list of all valid folders with ID:
sizeA = size(A,1);
for k=1:sizeA
    display(['Folder: ',A(k).name,' = No. ',num2str(k,'%02.0f')]);
end

toDoList=input(['There are: ',num2str(sizeA),' folders. Which do you want to analyze [all]: ']);

if isempty(toDoList)
    toDoList=1:sizeA;
end

% Now that the image processing is done do the statistics:
for i=toDoList
    folderName = A(i,1).name;
    
    display(['Working on folder: ',folderName]);
    
    cd([pwd,filesep,folderName]);
    
    % execute the migration script
    cellMig_part_2_statistics(binPix, minTrackLength, timeWindow, dframes, dt, showTrackerMovie, showEdgeMovie, 1);
    
    % go back to the super folder
    cd ..
    close all;
end