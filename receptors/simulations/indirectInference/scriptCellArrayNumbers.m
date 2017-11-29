currDir='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170501';

% cellInfoAllProbes=cell(100,5);
% count to know which row will be filled
temp = load(fullfile(currDir,filesep,'cellInfoAllSimProbes.mat'));
 cellInfoAllProbes=temp.cellInfoAllProbes;
% indexCount=1+size(cellInfoAllProbes,1);

%fill the columns with the info
dr=0.25:0.25:2;
for i=1:size(cellInfoAllProbes,1)

%second column is the receptor density
cellInfoAllProbes{i,6}=2:2:16;

%third column is the association probability

cellInfoAllProbes{i,7}=0.2:0.1:0.8;

%fourth column is the value of dissociation rate
cellInfoAllProbes{i,8}=dr(i);

%fifth column is the value of labeled fraction
cellInfoAllProbes{i,9}=0.1:0.1:0.6;

%load the cell array with the past info

% temp = load(fullfile(currDir,filesep,'cellInfoAllProbes.mat'));
% cellInfoAllProbesOld=temp.cellInfoAllProbes;

end

%add the new info in the Old file



%save updated cell array

save([currDir,filesep,'cellInfoAllProbes'],'cellInfoAllProbes','-v7.3');