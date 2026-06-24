function transplantMDProcesses(MD_new, MD_old, absBase)
% transplantMDProcesses
%
% Transplants processes and packages from MD_old into MD_new,
% then fixes all relative paths to absolute using absBase.
%
% USAGE:
%   transplantMDProcesses(MD_new, MD_old, absBase)
%
% INPUTS:
%   MD_new  : new MovieData with correct absolute paths (target)
%   MD_old  : old/backup MovieData with processes/packages (source)
%   absBase : absolute root to prepend to relative paths
%             e.g. '/mnt/nas/Collaborations/Celine_San_Diego/20260411'
%
% EXAMPLE:
%   transplantMDProcesses(MD_new, MD_old, ...
%       '/mnt/nas/Collaborations/Celine_San_Diego/20260411');

% ================================================================
% Step 1: Remove existing processes/packages from MD_new
% ================================================================
fprintf('[transplant] Clearing existing processes/packages from MD_new...\n');

for i = numel(MD_new.processes_):-1:1
    if ~isempty(MD_new.processes_{i})
        try
            MD_new.deleteProcess(MD_new.processes_{i});
        catch ME
            warning(ME.identifier, 'Could not delete process %d: %s', i, ME.message);
        end
    end
end

for i = numel(MD_new.packages_):-1:1
    if ~isempty(MD_new.packages_{i})
        try
            MD_new.deletePackage(MD_new.packages_{i});
        catch ME
            warning(ME.identifier, 'Could not delete package %d: %s', i, ME.message);
        end
    end
end

% ================================================================
% Step 2: Transplant processes from MD_old -> MD_new
% ================================================================
fprintf('[transplant] Transplanting %d processes...\n', numel(MD_old.processes_));

for i = 1:numel(MD_old.processes_)
    proc = MD_old.processes_{i};
    if isempty(proc), continue; end
    try
        MD_new.addProcess(proc);
        fprintf('  Added proc %d: %s\n', i, class(proc));
    catch ME
        warning(ME.identifier, 'Could not add process %d (%s): %s', i, class(proc), ME.message);
    end
end

% ================================================================
% Step 3: Transplant packages from MD_old -> MD_new
% ================================================================
fprintf('[transplant] Transplanting %d packages...\n', numel(MD_old.packages_));

for i = 1:numel(MD_old.packages_)
    pkg = MD_old.packages_{i};
    if isempty(pkg), continue; end
    try
        MD_new.addPackage(pkg);
        fprintf('  Added package %d: %s\n', i, class(pkg));
    catch ME
        warning(ME.identifier, 'Could not add package %d (%s): %s', i, class(pkg), ME.message);
    end
end

% ================================================================
% Step 4: Fix relative paths -> absolute in all processes
% ================================================================
fprintf('[transplant] Fixing relative paths (absBase: %s)...\n', absBase);

for i = 1:numel(MD_new.processes_)
    proc = MD_new.processes_{i};
    if isempty(proc), continue; end

    % Fix outFilePaths_
    try
        fp = proc.outFilePaths_;
        changed = false;
        for r = 1:size(fp,1)
            for c = 1:size(fp,2)
                if ischar(fp{r,c}) && ~isempty(fp{r,c}) && ~isAbsPath_(fp{r,c})
                    fp{r,c} = fullfile(absBase, fp{r,c});
                    changed = true;
                end
            end
        end
        if changed
            proc.setOutFilePaths(fp);
            fprintf('  Proc %d (%s): outFilePaths_ fixed\n', i, class(proc));
        end
    catch ME
        warning(ME.identifier, 'Proc %d outFilePaths_: %s', i, ME.message);
    end

    % Fix funParams_.OutputDirectory
    try
        p = proc.funParams_;
        if isfield(p, 'OutputDirectory') && ~isempty(p.OutputDirectory) ...
                && ~isAbsPath_(p.OutputDirectory)
            p.OutputDirectory = fullfile(absBase, p.OutputDirectory);
            proc.setPara(p);
            fprintf('  Proc %d (%s): funParams_.OutputDirectory fixed\n', i, class(proc));
        end
    catch ME
        warning(ME.identifier, 'Proc %d funParams_: %s', i, ME.message);
    end
end

% ================================================================
% Step 5: Save and verify
% ================================================================
MD_new.save();
fprintf('\n[transplant] Done. MD saved.\n');
fprintf('  outputDirectory_ : %s\n', MD_new.outputDirectory_);
fprintf('  processes_       : %d\n', numel(MD_new.processes_));
fprintf('  packages_        : %d\n', numel(MD_new.packages_));
fprintf('\n  Process summary:\n');
for i = 1:numel(MD_new.processes_)
    proc = MD_new.processes_{i};
    if isempty(proc)
        fprintf('    Proc %d: [empty]\n', i);
        continue;
    end
    % Check first output file exists
    fp = proc.outFilePaths_;
    firstFile = '';
    for r = 1:size(fp,1)
        for c = 1:size(fp,2)
            if ischar(fp{r,c}) && ~isempty(fp{r,c})
                firstFile = fp{r,c};
                break;
            end
        end
        if ~isempty(firstFile), break; end
    end
    fileOK = ~isempty(firstFile) && exist(firstFile, 'file') == 2;
    fprintf('    Proc %d: %-45s | file_ok=%d\n', i, class(proc), fileOK);
end
end

% ------------------------------------------------------------------
function tf = isAbsPath_(p)
tf = ~isempty(p) && (p(1) == '/' || (numel(p) > 1 && p(2) == ':'));
end