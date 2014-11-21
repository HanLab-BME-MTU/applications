function [exps] = whArrangeData(mainDirname)
% addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/LCCBcommon'));
% addpath(genpath('/work/gdanuser/azaritsky/TAU/UTSW/code'));
% addpath(genpath('/project/cellbiology/gdanuser/collab/assaf/TAU/UTSW/code/MET/woundHealing/utils'));
% addpath(genpath('/project/cellbiology/gdanuser/collab/assaf/TAU/UTSW/code/algs/'));
addpath(genpath('/project/cellbiology/gdanuser/collab/assaf/TAU/UTSW/code/utils/'));
addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/LCCBextern/bioformats/'));
addpath(genpath('/project/cellbiology/gdanuser/collab/assaf/TAU/UTSW/code/MET/woundHealing/utils/'));
% addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/Hall'));
% addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/utils'));
% addpath(genpath('/work/gdanuser/azaritsky/UTSW/code/algs'));
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
        
        if ~exist(dirname,'dir')
            unix(sprintf('mkdir %s',dirname));
        end
        
        imagesdir = [dirname '/images/'];
        
        if ~exist(imagesdir,'dir') && ~exist([dirname '/MF/'],'dir')
            unix(sprintf('mkdir %s',imagesdir));
            if (strcmp(ext, '.tif'))
                fname = [mainDirname name '.tif'];
                info = imfinfo(fname);
                nFrames = numel(info);
                for t = 1 : nFrames
                    I = imread(fname,t);
                    if size(I,3) > 1
                        I = I(:,:,1);
                    end
                    eval(['imwrite(I,''' [imagesdir pad(t,3) '.tif'''] ',''tif'')']);
                end
            else if (strcmp(ext, '.zvi'))                    
                    fname = [mainDirname name '.zvi'];
                    data = bfopen(fname);
                    images = data{1};
                    nFrames = size(images,1);
                    for t = 1 : nFrames
                        I = images(t,1);
                        I = I{:};
                        eval(['imwrite(I,''' [imagesdir pad(t,3) '.tif'''] ',''tif'')']);
                    end
                else if (strcmp(ext,'.lsm'))
                        fname = [mainDirname name '.lsm'];
                        stack = tiffread29(fname);
                        nFrames = length(stack);
                        for t = 1 : nFrames
                            data = stack(t).data;
                            if length(data) == 2
                                I = data{2};
                            else
                                I = data;
                            end
                            eval(['imwrite(I,''' [imagesdir pad(t,3) '.tif'''] ',''tif'')']);
                        end
                    end
                end
            end 
        else
            nFrames = length(dir(imagesdir)) - 2;
        end
        
        exps{iexps}.nFrames = nFrames;
        iexps = iexps + 1;
    end
end

end