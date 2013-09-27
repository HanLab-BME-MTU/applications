%sortTiffStacks() 

% Francois Aguet, 09/27/2013

function sortTiffStacks(varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('path', [], @(x) ischar(x) || isempty(x) || iscell(x));
ip.addParamValue('ChannelSpecifier', []);
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
    cellfun(@(x)sortTiffStacks(x, 'ChannelSpecifier', chSpec), stkpath);
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

for k = 1:N
    [~,dirName] = fileparts(stkList{k});
    [~,~] = mkdir([stkpath dirName]);
    
    movefile([stkpath stkList{k}], [stkpath dirName]);
end
