function [ output_args ] = GCAVisualsMontagingMovie(movieData,inputFolders,outputFolder,movieName)
%% GCAVisualsMontagingMovie: Generic file for montaging overlays 
% movieData : the movieData file where the overlays have been run
% previously
% inputFolders  : cell array of folders to montage 
% outputFolder : where the montage will be stored. 
% nameOutputMovie : movieName
% 
 if nargin<4 || isempty(movieName)
        movieName = 'movie';
 end
 
for fi = 1:movieData.nFrames_
    
    fstr = num2str(fi, '%.3d');
    
    % make a string
    forMont= cellfun(@(x) [' ' [x filesep fstr '.png'] ' '],inputFolders,...
        'uniformoutput',0); 
    filenamesStr = horzcat(forMont{:}); 
   
   
    outC = [outputFolder filesep 'frame' fstr '.png'];
    
    
    %Command for ImageMagick
    cmd = [' montage -tile x1 -geometry +5+5+0+0 -background "rgb(255,255,255)"' filenamesStr ...
        '  -compress lzw ' outC];
    system(cmd);
    %% TEST IF ON WINDOWS
    OSc = getenv('OS');
    if ~isempty(regexp(OSc,'Windows'),'ONCE')
        % if on a windows  machine test (if on linux do not need to crop -
        % windows just sucks)
        imgTest = imread(outC);
        [ny,nx] = size(imgTest);
        if (rem(ny,2)~=0)
            Y = '1' ;
        else
            Y = '0' ;
        end
        
        if rem(nx,2)~=0
            X = '1';
        else
            X = '0';
        end
        %
        cmd  = [' convert ' outC ' -crop ' num2str(nx) 'x' num2str(ny) '+' X '+' Y ' ' outC];
        % cmd = [' convert ' fname ' -shave ' toShaveX 'x' toShaveY ' -depth 8 ' fname];
        system(cmd);
        %% TAKE OUT BEFORE RELEASE - Line for personal use for some reason I was
        % having difficulty making ffmpeg run on windows without cd ing to the
        % original ffmpeg folder.
        cd(['C:\Users\Maria\Desktop\' ...
            'ffmpeg-20141210-git-ae81680-win64-static\ffmpeg-20141210-git-ae81680-win64-static\bin']);
    end



% name = projList{iProj,2};
%ffmpeg -i pathToFrames/frame_%0Nd.tif -r 15 -b 20000k filename.mp4
execute = ['ffmpeg -r 5 -i ' outputFolder filesep 'frame' '%03d.png' ...
    ' -crf 22 -pix_fmt yuv420p -b 20000k ' outputFolder filesep movieName '.mp4'];
system(execute)
cd(outputFolder)





end

