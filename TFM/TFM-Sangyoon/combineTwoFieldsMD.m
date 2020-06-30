function [] = combineTwoFieldsMD(MD,folderName1,folderName2)
%function [] = combineTwoFieldsMD(MD) uses MD to extract two displacement
%fields from folderName1 and folderName2, combine the two using
%combineTwoFields and store it as displacementField folder
%
% Sangyoon Han, June, 2020

% Load the two fields
if nargin<2
    folderName1='displacementField568';
    folderName2='displacementField640';
end

%main path
tPack = MD.getPackage(MD.getPackageIndex('TFMPackage'));
mainPath = tPack.outputDirectory_;
folderName1Full=[mainPath filesep folderName1];
folderName2Full=[mainPath filesep folderName2];
fileName = 'displField.mat';
displPath1 = [folderName1Full filesep fileName];
displPath2 = [folderName2Full filesep fileName];

% insepct if the folders are there
if exist(displPath1,'file')
    displObj = load(displPath1);
    displField1 = displObj.displField;
else
    folderName1Full=uigetdir(mainPath, 'Find the folder1');
    displPath1 = [folderName1Full filesep fileName];
    displObj = load(displPath1);
    displField1 = displObj.displField;
end

if exist(displPath2,'file')
    displObj = load(displPath2);
    displField2 = displObj.displField;
else
    folderName2Full=uigetdir(mainPath, 'Find the folder1');
    displPath2 = [folderName2Full filesep fileName];
    displObj = load(displPath2);
    displField2 = displObj.displField;
end

displField = combineTwoFields(displField1, displField2);

nFrames = MD.nFrames_;
[reg_grid,~,~,~]=createRegGridFromDisplField(displField,1.0,0);
displFieldShifted(nFrames)=struct('pos','','vec','');
for ii=1:nFrames
    % Shifted displField vector field
    [grid_mat,iu_mat, ~,~] = interp_vec2grid(displField(ii).pos, displField(ii).vec,[],reg_grid);
   
    [displFieldShiftedpos,displFieldShiftedvec, ~, ~] = interp_vec2grid(grid_mat+iu_mat, iu_mat,[],grid_mat); %1:cluster size
    pos = [reshape(displFieldShiftedpos(:,:,1),[],1) reshape(displFieldShiftedpos(:,:,2),[],1)]; %dense
    disp_vec = [reshape(displFieldShiftedvec(:,:,1),[],1) reshape(displFieldShiftedvec(:,:,2),[],1)]; 

    displFieldShifted(ii).pos = pos;
    displFieldShifted(ii).vec = disp_vec;
end

% Saving
displPack = tPack.processes_{2};
outputFile = displPack.outFilePaths_;
save(outputFile{1},'displField','displFieldShifted');
disp(['Done combining the two fields! It is overwritten in ' outputFile{1}])

