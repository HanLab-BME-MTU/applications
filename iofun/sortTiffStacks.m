%sortTiffStacks() 


% 1) Input directory contains TIFF stacks
%    Stacks are sorted into individual directories with the same file name.
% 2) List of TIFF stacks, including multiple channels. File names must contain a movie identifier.
%    Channels are sorted into their own directories. 
%    Example: cell1_1s_488.tif cell1_1s_561.tif cell2_1s_488.tif cell2_1s_561.tif,
% 3) 


% Francois Aguet, 09/27/2013

function sortTiffStacks(varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('path', []);%, @(x) ischar(x) || isempty(x) || iscell(x));
ip.addParamValue('ChannelNames', [], @iscell);
ip.addParamValue('MovieSelector', 'cell', @ischar);
ip.parse(varargin{:});
stkpath = ip.Results.path;
chSpec = ip.Results.ChannelNames;

if isempty(stkpath)
   stkpath = uigetdir('Select directory containing the STK files:'); 
   if (stkpath == 0)
       return;
   end
end

% Recursive call if input is cell array
if iscell(stkpath)
    cellfun(@(x)sortTiffStacks(x, chSpec), stkpath);
    return
end

stkpath = [stkpath filesep];
stkList = vdir(stkpath);

% if list of dirs: assume multi-channel, sort into sub-directories
if all([stkList.isdir])
    dirList = {stkList.name};
    for i = 1:numel(dirList)
        stkList = vdir([stkpath dirList{i}]);
        stkList = {stkList.name};
        stkList = stkList(~cellfun(@isempty, regexpi(stkList, '(\.tiff?|\.stk)$')));
        if ~isempty(chSpec)
            for c = 1:numel(stkList); % loop through channels
                [~,~,ext] = fileparts(stkList{c});
                destDir = [stkpath dirList{i} filesep chSpec{c} filesep];
                [~,~] = mkdir(destDir);
                movefile([stkpath dirList{i} filesep '*' chSpec{c} '*' ext], destDir);
            end
        else
            for c = 1:numel(stkList);
                [~,dname] = fileparts(stkList{c});
                destDir = [stkpath dirList{i} filesep dname filesep];
                [~,~] = mkdir(destDir);
                movefile([stkpath dirList{i} filesep stkList{c}], [destDir stkList{c}]);
            end
        end
        
    end
else % list of stacks
    stkList = {stkList(~[stkList.isdir]).name};
    idx = ~cellfun(@isempty, regexpi(stkList, '(\.tiff?|\.stk)$'));
    stkList = stkList(idx);
    
    N = length(stkList);
    if N==0
        fprintf('No TIFF files found in input directory.\n');
        return
    end
    
    
    % stack index/identifier following selector (i.e., 12 for 'Cell12_2s')
    movieID = str2double(regexpi(stkList, ['(?<=' ip.Results.MovieSelector ')\d+'], 'match', 'once'));
    
    % frame rate (if in file name)
    framerate = str2double(regexpi(stkList, '(?<=_)\d+(?=s)', 'match', 'once'));
    [movieID,idx] = unique(movieID);
    framerate = framerate(idx);
    
    nCh = numel(chSpec);
    
    N = numel(movieID);
    for k = N:-1:1
        if ~isnan(framerate(k))
            dirName = [ip.Results.MovieSelector num2str(movieID(k)) '_' num2str(framerate(k)) 's'];
        else
            dirName = [ip.Results.MovieSelector num2str(movieID(k))];
        end
        [~,~] = mkdir([stkpath dirName]);
        
        if nCh>0
            for c = 1:nCh
                % create sub-directory for each channel
                [~,~] = mkdir([stkpath dirName filesep chSpec{c}]);
                % move matching files to sub-directory
                movefile([stkpath ip.Results.MovieSelector num2str(movieID(k)) '*' chSpec{c} '*.tif*'],...
                    [stkpath dirName filesep chSpec{c}]);
            end
        else
            movefile([stkpath ip.Results.MovieSelector num2str(movieID(k)) '*.tif*'], [stkpath dirName]);
        end
    end
end


function d = vdir(path)
d = dir(path);
d = d(cellfun(@(i) ~strcmpi(i(1), '.'), {d.name}));
