expsFname1 = '/project/cellbiology/gdanuser/melanomaModel/Analysis/PrimaryMelanoma/PhaseContrast2DExperiments20140430_MD.mat';
expsFname2 = '/project/cellbiology/gdanuser/melanomaModel/Analysis/PrimaryMelanoma/PhaseContrast2DExperiments20140430.mat';

dirname1 = '/project/cellbiology/gdanuser/melanomaModel/Analysis/PrimaryMelanoma/';
dirname2 = '/work/gdanuser/azaritsky/UTSW/Data/Erik/POC_1min/';

load(expsFname1);
exps1 = metaData.experiments.fnames;
n1 = length(exps1);
load(expsFname2);
exps2 = metaData.experiments.fnames;
n2 = length(exps2);
clear metData;

for i = 1 : n2
    fname2 = exps2{i};
    for j = 1 : n1
        fname1 = exps1{j};
        if strcmp(fname2(1:6),fname1(1:6))
            curDirname1 = [dirname1 fname1 '/'];
            if ~exist(curDirname1,'dir')
                unix(sprintf('mkdir %s',curDirname1));
            end
            curDirname2 = [dirname2 fname2(1:6) '/'];
            %% copy results
            %             unix(sprintf('cp -R %s %s',[curDirname2 'results'],curDirname1));
            resultsDirname1 = [curDirname1 'results/'];
            resultsDirname2 = [curDirname2 'results/'];
            unix(sprintf('mkdir %s',resultsDirname1));
            filenames = dir(resultsDirname2);
            
            unix(sprintf('rm -rf %s',[resultsDirname1 '/*']));
            
            nfiles = length(filenames);
            for f = 3 : nfiles                
                filename = filenames(f).name;                
                filename2 = filename;
                
                % VERY SPECIFIC!!!!
                if length(filename2) < 12
                    continue;
                end
                                
                filename(7) = upper(filename(7));
                filename(11) = upper(filename(11));
                prefix = filename(1:14);
                ind = find(ismember(filename,'_'));
                suffix = filename(ind:end);
                task = filename(ind-2:ind-1);
                
                unix(sprintf('cp %s %s',[resultsDirname2 filename2],[resultsDirname1 prefix '_s' task suffix]));                
            end                        
            
            continue;
            
            %% copy all analysis files
            for k = 1 : 20
                taskDirname1 = [curDirname1 [fname1 '_s' pad(k,2)]];
                taskDirname2 = [curDirname2 [fname2 pad(k,2)]];
                unix(sprintf('cp -R %s %s',[taskDirname2 '/* ',taskDirname1 '/.']));                
            end
        end
    end
end