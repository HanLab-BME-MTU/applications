function analysisStruct = makiMakeAnalysisPlatformIndependent(analysisStruct)
%MAKIMAKEANALYSISPLATFORMINDEPENDENT makes an analysis structure compatible with current platform
%
% SYNOPSIS: analysisStruct = makiMakeAnalysisPlatformIndependent(analysisStruct)
%
% INPUT analysisStruct: maki analysisStruct with analysisStruct path and
%                       movie paths
%
% OUTPUT same as input with paths made compatible with current platform 
%
% REMARKS
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: kjaqaman
% DATE: 09-Aug-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%convert analysisStruct path
analysisStruct.filePath = makiPathDef(analysisStruct.filePath);

%get number of movies in analysisStruct
numMovies = size(analysisStruct.movies,1);

%convert the movie paths
for iMovie = 1 : numMovies
    analysisStruct.movies{iMovie,2} = makiPathDef(analysisStruct.movies{iMovie,2});
end
