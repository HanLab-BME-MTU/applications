function [ output_args ] = GCAVisualsMontagingMovie(movieData,inputFolders,outputFolder,varargin)
%% GCAVisualsMontagingMovie: Generic file for montaging overlays 
% movieData : the movieData file where the overlays have been run
% previously
% inputFolders  : cell array of folders to montage 
% outputFolder : where the montage will be stored. 
% nameOutputMovie : movieName
% 

%%Input check
ip = inputParser;

ip.CaseSensitive = false;
ip.addRequired('inputFolders'); 
ip.addRequired('outputFolder'); 
ip.addParameter('movieName','movie');
ip.addParameter('makeMontage',true); 
ip.addParameter('runffmpeg',true); 
ip.addParameter('tiling', 'x1'); %  default single row, 2x2 for square tiling 'x1'

ip.parse(inputFolders,outputFolder,varargin{:}); 
 
movieName = ip.Results.movieName; 

if ip.Results.makeMontage
for fi = 1:movieData.nFrames_
    
    fstr = num2str(fi, '%.3d');
    
    % make a string
    forMont= cellfun(@(x) [' ' [x filesep fstr '.png'] ' '],inputFolders,...
        'uniformoutput',0); 
    filenamesStr = horzcat(forMont{:}); 
   
   
    outC = [outputFolder filesep 'frame' fstr '.png'];
    
    
    %Command for ImageMagick
%    cmd = [' montage -tile x1 -geometry +5+5+0+0 -background "rgb(255,255,255)"' filenamesStr ...
%     '  -compress lzw ' outC];

 cmd = [' montage -tile ' ip.Results.tiling ' -geometry +5+5+0+0 -background "rgb(255,255,255)"' filenamesStr ...
    '  -compress lzw ' outC];
%       cmd = [' montage -tile 2x2 -geometry +5+5+0+0 -background "rgb(255,255,255)"' filenamesStr ...
%     '  -compress lzw ' outC];
    
    
    system(cmd);
end 
end 
    %% TEST IF ON WINDOWS
    OSc = getenv('OS');
    if ~isempty(regexp(OSc,'Windows','ONCE'))
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

if ip.Results.runffmpeg

% name = projList{iProj,2};
%ffmpeg -i pathToFrames/frame_%0Nd.tif -r 15 -b 20000k filename.mp4
execute = ['ffmpeg -r 5 -i ' outputFolder filesep 'frame' '%03d.png' ...
    ' -crf 22 -pix_fmt yuv420p -b 20000k ' outputFolder filesep movieName '.mp4']; % note this codec I am having trouble with on my windows pc
system(execute)

execute = ['ffmpeg -r 5 -i ' outputFolder filesep 'frame' '%03d.png' ...
    ' -b 2000k ' outputFolder filesep movieName  '.wmv'];
system(execute);

execute = ['ffmpeg -r 5 -i ' outputFolder filesep 'frame' '%03d.png' ...
    ' -crf 22 -b 20000k' outputFolder filesep movieName '2.mp4'];
system(execute);


cd(outputFolder)

end 




end

