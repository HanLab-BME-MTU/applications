A = dir();
sizeA = size(A,1);

% if the folder has files that need to be sorted into folders first, use
% the following:

positions = 12:23;

for p = positions
    newFolder = sprintf('%s%.2d','xy',p);
    mkdir(newFolder)
    
    %     for s=1:sizeA
    %         fileName = A(s,1).name;
    source = sprintf('%s%.2d%s','*xy',p,'*.tif');
    movefile(source,newFolder)
    %     end
end

A = dir();
sizeA = size(A,1);

% these are the Folder Names (FN):
%FN_data    = 'data';
FN_Phase   = 'Phase';  % w1
FN_Nuclei   = 'Nuclei';  % w2

for i=1:sizeA
    folderName = A(i,1).name
    parentFolder = sprintf('%s%s%s',pwd,filesep,folderName,filesep);
    
%     mkdir(parentFolder,FN_data);
%     subFolder = sprintf('%s%s',parentFolder,filesep,FN_data,filesep);    
    
    % treat the files with pattern *c1*:
    sourceFilesCell = sprintf('%s%s',parentFolder,'*c1*');
    if ~isempty(dir(sourceFilesCell))
        % mkdir(subFolder, FN_Phase);
        mkdir(parentFolder, FN_Phase);
        %destFilesCell = sprintf('%s%s',parentFolder,FN_data,filesep,FN_Phase);
        destFilesCell = sprintf('%s%s',parentFolder,FN_Phase);
        movefile(sourceFilesCell,destFilesCell);
    else
        display('Pattern: *w1* was not found')
    end
    
    % treat the files with pattern *w2*:
    sourceFilesNuclei = sprintf('%s%s',parentFolder,'*c2*');
    if ~isempty(dir(sourceFilesNuclei))
        %mkdir(subFolder, FN_Nuclei);
        mkdir(parentFolder, FN_Nuclei);
        %destFilesNuclei = sprintf('%s%s',parentFolder,FN_data,filesep,FN_Nuclei);
        destFilesNuclei = sprintf('%s%s',parentFolder,FN_Nuclei);
        movefile(sourceFilesNuclei,destFilesNuclei);
    else
        display('Pattern: *w2* was not found')
    end
    
 
end
%%
% for i=5:sizeA
%     folderName = A(i,1).name;
%     parentFolder = sprintf('%s%s%s','/home/mrn8/orchestra/groups/lccb-mechanics/analysis/2010_05_11_TFM_MDA231_tobeanalyzed/',folderName,'/');
%     mkdir(parentFolder,'data');
%     
%     subFolder = sprintf('%s%s',parentFolder,'/data/');
%     mkdir(subFolder,'Beads');
%     mkdir(subFolder,'Cell');
%     mkdir(subFolder,'Reference Frame');
% %     mkdir(subFolder,'Ecad'); % Comment as needed
%     mkdir(subFolder,'CAAX'); % Comment as needed
%     
%     sourceFilesRef = sprintf('%s%s',parentFolder,'*REF*'); % Important to do this first since the filenames also contain w1DIC and w2647 which we use later
%     destFilesRef = sprintf('%s%s',parentFolder,'data/Reference Frame/');
%     movefile(sourceFilesRef,destFilesRef);
%     
%     sourceFilesMovies = sprintf('%s%s',parentFolder,'*.mp4'); % Important to do this first since the filenames also contain w1DIC and w2647 which we use later
%     destFilesMovies = sprintf('%s%s',parentFolder,'data/');
%     movefile(sourceFilesMovies,destFilesMovies);
%     
%     sourceFilesCell = sprintf('%s%s',parentFolder,'*w1*');
%     destFilesCell = sprintf('%s%s',parentFolder,'data/Cell/');
%     movefile(sourceFilesCell,destFilesCell);
%     
%     sourceFilesBeads = sprintf('%s%s',parentFolder,'*w2*');
%     destFilesBeads = sprintf('%s%s',parentFolder,'data/Beads/');
%     movefile(sourceFilesBeads,destFilesBeads);
%     
% %     %The following is as needed
% %     sourceFilesEcad = sprintf('%s%s',parentFolder,'*w3*');
% %     destFilesEcad = sprintf('%s%s',parentFolder,'data/Ecad/');
% %     movefile(sourceFilesEcad,destFilesEcad);
%     
%      %The following is as needed
%     sourceFilesEcad = sprintf('%s%s',parentFolder,'*w4*');
%     destFilesEcad = sprintf('%s%s',parentFolder,'data/CAAX/');
%     movefile(sourceFilesEcad,destFilesEcad);
% end