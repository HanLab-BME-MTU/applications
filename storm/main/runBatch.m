function runBatch(nCores)

% Define the root path
queuePath = [getStormPath() filesep '_queue' filesep];

% Host list
hostlist = ['clarinet001-032.orchestra ' ...
'clarinet001-033.orchestra ' ... % 8-cores
'clarinet001-034.orchestra ' ...
'clarinet001-035.orchestra ' ...
'clarinet001-036.orchestra ' ...
'clarinet001-037.orchestra ' ...
'clarinet001-038.orchestra ' ...
'clarinet001-039.orchestra ' ...
...
'clarinet002-047.orchestra ' ... % 12-cores
'clarinet002-048.orchestra ' ...
...% 'clarinet002-049.orchestra ' ...
...% 'clarinet002-050.orchestra ' ...
...% 'clarinet002-051.orchestra ' ...
...% 'clarinet002-052.orchestra ' ...
'clarinet002-053.orchestra'];

% Parameters
if ~nargin
    nCoresMin = 8;
    nCoresMax = 12;
else
    nCoresMin = nCores;
    nCoresMax = nCores;
end

if nCoresMin < 10 && nCoresMax < 10
    fprintf('######################## BATCH %.0f-%.0f #####################\n',nCoresMin,nCoresMax);
elseif nCoresMin >= 10 && nCoresMax >= 10
    fprintf('######################## BATCH %.0f-%.0f ###################\n',nCoresMin,nCoresMax);
else
    fprintf('######################## BATCH %.0f-%.0f ####################\n',nCoresMin,nCoresMax);
end

% List input directory
listing = dir(queuePath);
queueIdx = 0;
dirName = cell(numel(listing)-2,1);
for i=3:numel(listing)
    % Find all the directories with an appended dash
    if strcmp(listing(i).name(end),'-')
        queueIdx = queueIdx + 1;
        
        src = [queuePath listing(i).name];
        newName = [listing(i).name(1:end-1) '='];
        dest = [queuePath newName];
        movefile(src,dest,'f');
        
        dirName{queueIdx} = newName;
    end
end
dirName = dirName(1:queueIdx);
nDirName = numel(dirName);

if ~isempty(dirName)
    % Set up the scheduler
    jm = findResource('scheduler','type','lsf');
    set(jm,'ClusterMatlabRoot','/opt/matlab');
    set(jm, 'DataLocation','~/jobs');
    set(jm, 'SubmitArguments', ['-n ' num2str(nCoresMin) ' -R "rusage[mem=4000:matlab_dc_lic=1] && span[hosts=1]" -q danuser_7d']);
    % set(jm, 'SubmitArguments', ['-n ' num2str(nCoresMin) ' -R "rusage[mem=4000:matlab_dc_lic=1] && span[hosts=1]" -q danuser_7d -m ' hostlist]);
    
    % Process all the inputs
    job = cell(numel(dirName),1);
    
    for file=1:numel(dirName)
        dirPath = [queuePath dirName{file} filesep];
        
        % Create the job
        job{file} = createMatlabPoolJob(jm);
        set(job{file},'MaximumNumberOfWorkers',nCoresMax);
        set(job{file},'MinimumNumberOfWorkers',nCoresMin);
        
        % Create the task
        createTask(job{file},@algorithm,1,{dirPath});
        
        % Submit the job
        submit(job{file});
        
        % Display status
        fprintf('Main: Submitting job %3u/%u\n',file,nDirName);
    end
    
    % 0:no 1:pending 2:queued 3:running 4:finished 5:failed 6:destroyed 7:unavailable
    
    oldStates = zeros(numel(dirName),1);
    
    while ~all(oldStates > 3)
        pause(5);
        
        for file=1:numel(dirName)
            stateString = get(job{file},'State');
            if strcmp(stateString,'pending')
                stateIdx = 1;
            elseif strcmp(stateString,'queued')
                stateIdx = 2;
            elseif strcmp(stateString,'running')
                stateIdx = 3;
            elseif strcmp(stateString,'finished')
                stateIdx = 4;
            elseif strcmp(stateString,'failed')
                stateIdx = 5;
            elseif strcmp(stateString,'destroyed')
                stateIdx = 6;
            elseif strcmp(stateString,'unavailable')
                stateIdx = 7;
            else
                stateIdx = 0;
            end
            
            if stateIdx > oldStates(file)
                if strcmp(stateString,'finished')
					exitflag = getAllOutputArguments(job{file});
                    if numel(exitflag)
                        stateString = sprintf('FINISHED %s',exitflag{1});
                    else
                        stateString = sprintf('FINISHED Crashed?!?!?');
                    end
                elseif strcmp(stateString,'failed')
                    stateString = 'FAILED';
                elseif strcmp(stateString,'destroyed')
                    stateString = 'DESTROYED';
                end
                fprintf('%d: %s %s %s\n',file,dirName{file},datestr(now,'HH-MM-SS'),stateString);
                oldStates(file) = stateIdx;
            end
        end
    end
            
else
    disp('Main: Nothing to process!');
end

disp('######################## DONE ##########################');

end