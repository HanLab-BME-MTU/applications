%sortTiffStacks() 

% Francois Aguet, 09/27/2013

function sortTiffStacks(varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('path', [], @(x) ischar(x) || isempty(x) || iscell(x));
ip.addOptional('ChannelSpecifier', [], @iscell);
ip.addParamValue('MovieSelector', 'cell', @ischar);
ip.parse(varargin{:});
stkpath = ip.Results.path;
chSpec = ip.Results.ChannelSpecifier;

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
stkList = dir(stkpath);
stkList = {stkList(~[stkList.isdir]).name};
idx = ~cellfun(@isempty, regexpi(stkList, '(\.tiff?|\.stk)$'));
stkList = stkList(idx);

N = length(stkList);
if N==0
    fprintf('No TIFF files found in input directory.\n');
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
            movefile([stkpath ip.Results.MovieSelector num2str(movieID(k)) '*' chSpec{c} '*.*'],...
                [stkpath dirName filesep chSpec{c}]);
        end
    else
        movefile([stkpath ip.Results.MovieSelector num2str(movieID(k)) '*.*'], [stkpath dirName]);
    end
end
