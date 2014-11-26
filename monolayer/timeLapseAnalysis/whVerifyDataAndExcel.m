function [] = whVerifyDataAndExcel()
addpath(genpath('/work/gdanuser/azaritsky/TAU/UTSW/code'));
addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/Hall'));
addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/utils'));
addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/algs'));

mainDirname = '/work/gdanuser/azaritsky/UTSW/Data/Hall/Data20140224/';
expsData = getExpsData(mainDirname);
load([mainDirname 'GTPasesScreenMetaData20140121_new.mat']); % exps
expsExcel = exps;
clear exps;

hits = false(1,length(expsData));
for iexcel = 1 :length(expsExcel.fnames)
    curExcel = expsExcel.fnames{iexcel};
    found = false;
    for iexp = 1 : length(expsData)
        curExp = expsData{iexp}.name;
        
        if strcmp(curExcel,curExp)  
            if hits(iexp)
                % should not happend
                fprintf(sprintf('%s twice in experiment files\n',curExp));
                break;
            end
            hits(iexp) = true;
            found = true;
            break;
        end
    end
    if ~found
        fprintf(sprintf('%s not found in experiment files\n',curExcel));
    end
end

inds = find(~hits);
for iexp = 1 : length(inds)
    fprintf(sprintf('%s no match in excel!\n',expsData{iexp}.name));
end

end


%%
function [exps] = getExpsData(mainDirname)

filenames = dir(mainDirname);
nfiles = length(filenames);

iexps = 1;

for i = 1 : nfiles
    filename = filenames(i).name;
    
    [pathstr, name, ext] = fileparts(filename);
    
    if (strcmp(ext, '.tif') || strcmp(ext, '.zvi')|| strcmp(ext, '.lsm'))        
        dirname = [mainDirname name];
        
        exps{iexps}.name = name;
        exps{iexps}.ext = ext;
                
        iexps = iexps + 1;
    end
end

end