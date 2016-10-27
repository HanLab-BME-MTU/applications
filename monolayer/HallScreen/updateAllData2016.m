% Arrange raw data for 20160314
function [] = updateAllData2016()


workdir = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/';

allDataDir = [workdir 'allData/'];

dataDirs = {[workdir '/Data20150928/'],[workdir '/Data20160201/'],[workdir '/Data20160314/']};

nDataDirs = length(dataDirs);

for i = 1 : nDataDirs
    curdir = dataDirs{i};
    unix(sprintf('cp -R %s/* %s/.',curdir,allDataDir)); 
    fprintf(sprintf('copied:\n%s\nto\n%s\n\n',curdir,allDataDir));
end

end