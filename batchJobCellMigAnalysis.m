% browse into the folder that contains all the subfolders and execute this
% script.

% set the parameters for the data analysis:
display('Now set the parameters for the statistics:');
r   =input('Specify the minimal nuclei radius    r=[ 5pix]: ');
rMin=input('Specify the minimal search radius rMin=[    r]: ');
rMax=input('Specify the maximal search radius rMax=[50pix]: ');
smoothFac=input('Specify smoothing factor for smoother edges [2]: ');

if isempty(r)
    r=5;
end
if isempty(rMin)
    rMin=r;
end
if isempty(rMax)
    rMax=50;
end
if isempty(smoothFac)
    smoothFac=2;
end

resultDir=pwd;
if ~isdir(resultDir)
    mkdir(resultDir);
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

for i=toDoList
    folderName = A(i,1).name;
    
    display(['Working on folder: ',folderName]);
    
    cd([pwd,filesep,folderName]);
    
    % execute the migration script
    cellMig_part_1_analysis(r,rMin,rMax,smoothFac,[],[],0);
    
    % go back to the super folder
    cd ..
    close all;
end