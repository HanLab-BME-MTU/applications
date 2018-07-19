function [] = pcListTasksAll(dataDir)
% Outputs list of folders to a text file at MetaData folder
% Assaf Zaritsky, October 2017

close all;

if nargin  == 0
    dataDir = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/All';
end

assert(logical(exist(dataDir,'dir')));

outFname = [dataDir filesep '..' filesep 'MetaData' filesep 'allTaskFolders.txt'];
fid = fopen(outFname,'w');

filenamesExp = dir(dataDir);
nfilesExp = length(filenamesExp);

for iexp = 3 : nfilesExp
    filenameExp = filenamesExp(iexp).name;
    
    [pathstr, nameExp, ext] = fileparts(filenameExp);
    expDir = [dataDir filesep nameExp];
    if exist(expDir,'dir')
        
        filenamesTask = dir(expDir);
        nfilesTask = length(filenamesTask);
        
        for itask = 3 : nfilesTask
            filenameTask = filenamesTask(itask).name;
            
            [pathstr, nameTask, ext] = fileparts(filenameTask);
            taskDir = [expDir filesep nameTask];
            if exist(taskDir,'dir')
                fprintf(fid,sprintf('%s\n',taskDir));
            end
        end
    end    
end
fclose(fid);
end
