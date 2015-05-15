function [ output_args ] = gcaVisualsCropImageToMakeMovieOnWindows(movieData,visType)
%gcaVisualsCropImageToMakeMovieOnWindows. 
% for some reason windows is a bit picky about the dimensions of your image
% when running through ffmpeg (at least with my set up) therefore here we 
% simply shave to make the proportions of each image such that windows is happy. 
% we can then make the movies in ffmpeg no problem 
% 
%% START 
[ny,nx] = movieData.imSize_; 
cropDir = [movieData.outputDirectory_ filesep 'Visualization_Overlays' ... 
    filesep visType]; % for now juat assume you are cropping 
if ~isdir(cropDir) 
    mkdir(cropDir) 
end 

for fi = 1:movieData.nFrames_
    fstr = num2str(fi, '%.3d');
      
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
    nameC = [cropDir filesep num2str(fstr) '.png']; 
    cmd  = [' convert ' outC ' -crop ' num2str(nx) 'x' num2str(ny) '+' X '+' Y ' ' outC];
    % cmd = [' convert ' fname ' -shave ' toShaveX 'x' toShaveY ' -depth 8 ' fname];
    system(cmd);
end

