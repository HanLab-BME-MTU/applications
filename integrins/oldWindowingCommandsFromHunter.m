
%setup project
movieData = setupMovieData;

%view images
imageViewer(movieData);

%create masks - compare to PANDA
%can be inserted manually if PANDA is better
movieData = thresholdMovie(movieData,'MaxJump',0.2,'GaussFilterSigma',1);
movieData = refineMovieMasks(movieData,'ClosureRadius',[0 5]);

%makes movie with masks
makeMaskMovie(movieData);

%make windows
%these old windows always have straight edges between the contours
%to feed them into inpolygon, the innerBorder must be given in reverse
movieData = getMovieContours(movieData,(0:10:1000),[],1,[],[],2);
movieData = getMovieWindows(movieData,'e',10,[],1,10);

%create movie with windows
%if protrusion vectors are there, they are displayes as well
makeWindowTestingMovie(movieData);

%get protrusion vectors - compare to PANDA
movieData = getMovieProtrusion(movieData,1);

%to insert protrusion vectors manually if PANDA is better ...
movieData.protrusion.directory = pwd;
movieData.protrusion.fileName = 'protrusion.mat';
movieData.protrusion.nfileName = 'normal_matrix.mat';

%to average protrusion vectors in each window
movieData = getMovieProtrusionSamples(movieData);

%display data as image - this is a matlab function
imagesc(protrusionSamples.averageNormalComponent);

%center image color-coding
caxis([-2 2])
 
%average intensity in each window
%use this function as a template to get my integrin parameters in each
%window
movieData = getMovieActivitySamples(movieData);

%% Final sequence of commands:

movieData = setupMovieData;
movieData = thresholdMovie(movieData,'MaxJump',0.2,'GaussFilterSigma',1);
movieData = refineMovieMasks(movieData,'ClosureRadius',[0 5]);
makeMaskMovie(movieData);

movieData = getMovieProtrusion(movieData,1);
movieData = getMovieContours(movieData,(0:5:1000),[],1,[],[],2);
movieData = getMovieWindows(movieData,'e',5,[],1,2);
movieData = getMovieProtrusionSamples(movieData);
close all
makeWindowTestingMovie(movieData);
