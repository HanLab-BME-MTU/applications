A = dir('well503_*');
sizeA = size(A,1);

% these are the Folder Names (FN):
FN_data    = 'data';
FN_Phase   = 'Cells';  % w1
FN_Beads   = 'Beads';  % w2
FN_Xtra1st = 'Ecad';  % w3 This is actin or myosin or Ecad?!
FN_Xtra2nd = 'Myos';  % w4 This is actin or myosin or Ecad?!
FN_RefFrame= 'Reference Frame';

for i=1:sizeA
    folderName = A(i,1).name
    parentFolder = sprintf('%s%s%s',pwd,filesep,folderName,filesep);
    
    mkdir(parentFolder,FN_data);
    subFolder = sprintf('%s%s',parentFolder,filesep,FN_data,filesep);    
    
    % treat the reference frames:
    sourceFilesRef = sprintf('%s%s',parentFolder,'*REF*'); % Important to do this first since the filenames also contain w1DIC and w2647 which we use later
    if ~isempty(dir(sourceFilesRef))
        mkdir(subFolder, FN_RefFrame);
        destFilesRef = sprintf('%s%s',parentFolder,FN_data,filesep,FN_RefFrame);
        movefile(sourceFilesRef,destFilesRef);
    else
        display('No reference frames found!')
        display('Go and get them, nothing has been done!')
        return;
    end    
    
    % treat the movie files:
    sourceFilesMovies = sprintf('%s%s',parentFolder,'*.mp4'); % Important to do this first since the filenames also contain w1DIC and w2647 which we use later
    if ~isempty(dir(sourceFilesMovies))
        destFilesMovies = sprintf('%s%s',parentFolder,FN_data,filesep);
        movefile(sourceFilesMovies,destFilesMovies);
    else
        display('No .mp4 movies found!')
    end    
    
    % treat the files with pattern *w1*:
    sourceFilesCell = sprintf('%s%s',parentFolder,'*w1*');
    if ~isempty(dir(sourceFilesCell))
        mkdir(subFolder, FN_Phase);
        destFilesCell = sprintf('%s%s',parentFolder,FN_data,filesep,FN_Phase);
        movefile(sourceFilesCell,destFilesCell);
    else
        display('Pattern: *w1* was not found')
    end
    
    % treat the files with pattern *w2*:
    sourceFilesBeads = sprintf('%s%s',parentFolder,'*w2*');
    if ~isempty(dir(sourceFilesBeads))
        mkdir(subFolder, FN_Beads);
        destFilesBeads = sprintf('%s%s',parentFolder,FN_data,filesep,FN_Beads);
        movefile(sourceFilesBeads,destFilesBeads);
    else
        display('Pattern: *w2* was not found')
    end
    
    % treat the files with pattern *w3*:
    sourceFilesXtra1st = sprintf('%s%s',parentFolder,'*w3*');
    if ~isempty(dir(sourceFilesXtra1st))
        mkdir(subFolder, FN_Xtra1st);
        destFilesEcad = sprintf('%s%s',parentFolder,FN_data,filesep,FN_Xtra1st);
        movefile(sourceFilesXtra1st,destFilesEcad);
    else
        display('Pattern: *w3* was not found')
    end
    
    % treat the files with pattern *w4*:
    sourceFilesXtra2nd = sprintf('%s%s',parentFolder,'*w4*');
    if ~isempty(dir(sourceFilesXtra2nd))
        mkdir(subFolder, FN_Xtra2nd);
        destFilesEcad = sprintf('%s%s',parentFolder,FN_data,filesep,FN_Xtra2nd);
        movefile(sourceFilesXtra2nd,destFilesEcad);
    else
        display('Pattern: *w4* was not found')
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