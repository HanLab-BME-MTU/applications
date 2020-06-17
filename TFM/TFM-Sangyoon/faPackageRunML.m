function [] = faPackageRunML(MLPath,processesToRun)
%function [] = faPackageRunML(MLPath,processesToRun) runs faPackageRun per
%movieData in ML.
% Sangyoon Han, May 2020

if isa(MLPath,'MovieList')
    ML=MovieList.load(MLPath.getFullPath); %,'askUser',false,'askUserChannel',false);
else
    ML=MovieList.load(MLPath,'askUser',false); %,'askUserChannel',false);
end

nMovies = numel(ML.movieDataFile_);
MDAll = ML.movies_;
for ii=1:nMovies
    curMD = MDAll{ii};
    if nargin<2
        faPackageRun(curMD)
    else
        faPackageRun(curMD, processesToRun)
    end
end
ML.save
end
