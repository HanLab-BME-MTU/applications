% Arrange raw data for 20160314
function [] = whMoveRawData20160314()

workdir = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Data20160314/';
transferDir = [workdir 'TransferYunYu/'];

% cd(transferDir);

filenames = dir(transferDir);
nfiles = length(filenames);

for i = 1 : nfiles
    filename = filenames(i).name;
    [pathstr, name, ext] = fileparts(filename);
    
    if (strcmp(ext,'.zip'))
        curSourceDir = [transferDir name filesep];
        curTargetDirPrefix = [workdir name];
        
        for s = 1 : 99
            tic;
            
            d = dir(sprintf('%s*s%02d*.tif',curSourceDir,s));
            if isempty(d) 
                break;
            end
            
            curdir = [curTargetDirPrefix sprintf('_s%02d',s) filesep];
            imgdir = [curdir 'images' filesep];
            if ~exist(curdir,'dir')
                unix(sprintf('mkdir %s',curdir));
                unix(sprintf('mkdir %s',imgdir));
            else
                warning('folder %s exists!',curdir);
            end
            unix(sprintf('cp %s*s%02d*.tif %s',curSourceDir,s,imgdir));
            
            imgnames = dir(imgdir);
            nimages = length(imgnames);
                        
            cd(imgdir);
            
            for im = 1 : nimages
                imgname = imgnames(im).name;    
                [impathstr, imname, imext] = fileparts(imgname);
                if (strcmp(imext,'.tif'))
                    %                     unix(sprintf('mv %s %s',[imgdir imname imext],[imgdir imname(end-2:end) imext]));
                    movefile([imname imext],[imname(end-2:end) imext]);
                end
            end
            tt = toc;
            fprintf(sprintf('%s: %d (sec) \n',[name sprintf('_s%02d',s)],tt));
        end
        fprintf(sprintf('\n\n%s: %d locations \n',name,s-1));
        tic
    end
end


% iexps = 1;
% 
% nNew = 0;
% 
% for i = 1 : nfiles
%     filename = filenames(i).name;
%     
%     [pathstr, name, ext] = fileparts(filename);
%     
%     if (strcmp(ext, '.tif') || strcmp(ext, '.zvi')|| strcmp(ext, '.lsm'))        
%         dirname = [mainDirname name];
%         
%         exps{iexps}.name = name;
%         exps{iexps}.ext = ext;
%         
%         if ~exist(dirname,'dir')
%             unix(sprintf('mkdir %s',dirname));
%         end
%         
%         imagesdir = [dirname '/images/'];
%         
%         if ~exist(imagesdir,'dir') && ~exist([dirname '/MF/'],'dir')
%             nNew = nNew + 1;
%             unix(sprintf('mkdir %s',imagesdir));
%             if (strcmp(ext, '.tif'))
%                 fname = [mainDirname name '.tif'];
%                 info = imfinfo(fname);
%                 nFrames = numel(info);
%                 for t = 1 : nFrames
%                     I = imread(fname,t);
%                     if size(I,3) > 1
%                         I = I(:,:,1);
%                     end
%                     eval(['imwrite(I,''' [imagesdir pad(t,3) '.tif'''] ',''tif'')']);
%                 end
%             else if (strcmp(ext, '.zvi'))                    
%                     fname = [mainDirname name '.zvi'];
%                     data = bfopen(fname);
%                     images = data{1};
%                     nFrames = size(images,1);
%                     for t = 1 : nFrames
%                         I = images(t,1);
%                         I = I{:};
%                         eval(['imwrite(I,''' [imagesdir pad(t,3) '.tif'''] ',''tif'')']);
%                     end
%                 else if (strcmp(ext,'.lsm'))
%                         fname = [mainDirname name '.lsm'];
%                         stack = tiffread29(fname);
%                         nFrames = length(stack);
%                         for t = 1 : nFrames
%                             data = stack(t).data;
%                             if length(data) == 2
%                                 I = data{2};
%                             else
%                                 I = data;
%                             end
%                             eval(['imwrite(I,''' [imagesdir pad(t,3) '.tif'''] ',''tif'')']);
%                         end
%                     end
%                 end
%             end 
%         else
%             nFrames = length(dir(imagesdir)) - 2;
%         end
%         
%         exps{iexps}.nFrames = nFrames;
%         iexps = iexps + 1;
%     end
% end
% 
% fprintf(sprintf('\n\n new files = %d\n\n',nNew));
% 
% end