% This function provides the ouput of the 'result' of cmeAnalysis as an
% excel file in the directory chosen by the user. 
% 1. Load 'result' or run 'cmeAnalysis'
% 2. Call the function: xxx = cmeOutput(result), xxx is the output name
% 3. Choose preferred output folder
% 4. If the output file already exists, select replace or not
%=============================================================
% Required Matlab version: 2017
%=============================================================
% Optionally, the name of the excel file can be manually changed
%=============================================================
% Author: Xinxin Wang, Danuser Lab

function pwd_file = cmeOutput(res, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('res', @(x) isstruct(x));
ip.addOptional('outputFilename', 'res', @ischar);
ip.addOptional('resultDir', [], @(x) ischar(x));
ip.parse(res, varargin{:});

filename = [ip.Results.outputFilename '.xlsx'];
resDir = ip.Results.resultDir;

if isempty(resDir)
    resDir = uigetdir(pwd, 'Select output folder:');
    if resDir==0
        return
    else
        resDir = [resDir filesep];
    end
end
if ~strcmp(resDir(end), filesep)
    resDir = [resDir filesep];
end
pwd_file = [resDir filename];
cd(resDir);

if exist(filename, 'file')
  disp('Warning: file already exists:');
  disp(filename);
  askReplaceOrNot = false;  
  while askReplaceOrNot == false
  prompt = 'Do you want to replace? Y/N: ';
  strReplaceOrNot = input(prompt,'s');
  if (strReplaceOrNot == 'N') || (strReplaceOrNot == 'n')
    return
  elseif (strReplaceOrNot == 'Y') || (strReplaceOrNot == 'y')
    delete(filename);
    askReplaceOrNot = true;   
  end
  
  end
  disp('Output file name:');
  disp(filename);
else
  disp('Output file name:');
  disp(filename);
end

meanLftHistCCP = res.lftRes.meanLftHistCCP';

Lft = 5:size(res.lftRes.meanLftHistCCP,2)+4;
Lft = Lft';
cumMeanLftHistCCP = cumsum(meanLftHistCCP);
T = table(Lft,meanLftHistCCP,cumMeanLftHistCCP);
writetable(T,filename,'Sheet',1);


num_movie = size(res.lftRes.cellArea,1);
cumLftHistCCP = cumsum(res.lftRes.lftHistCCP(1,:));
for i = 2: num_movie
    cumLftHistCCP = cat(1,cumLftHistCCP,cumsum(res.lftRes.lftHistCCP(i,:)));
end

LT_50 = zeros(1,num_movie);
for i = 1: num_movie
val_50 = 0.5;
tmp = abs(cumLftHistCCP(i,:)-val_50);
[idx idx] = min(tmp);
%LT_50(i) = cumLftHistCCP(i,idx)
LT_50(i) = Lft(idx);
end
LT_50_mean = mean(LT_50);
LT_50_std = std(LT_50);
T = table(LT_50_mean,LT_50_std);

writetable(T,filename,'Sheet',1, 'Range','E1');

cellNum = 1:size(res.lftRes.cellArea,1);
cellNum = cellNum';
cellArea = res.lftRes.cellArea;
%pctCCP = res.lftRes.pctCCP;
%pctCS = res.lftRes.pctCS;
%pctVisit = res.lftRes.pctVisit;
initDensityCCP = res.lftRes.initDensityCCP(:,1);
initDensityCCP_all = res.lftRes.initDensityIa(:,1);
initDensityCLS = initDensityCCP_all-initDensityCCP;
persistentDensity = res.lftRes.persistentDensity;
LT_50 = LT_50';
T = table(cellNum,cellArea,initDensityCCP,initDensityCLS,persistentDensity,LT_50);

writetable(T,filename,'Sheet',1, 'Range','H1');

