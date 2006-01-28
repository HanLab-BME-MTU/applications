function imarisHandle = imarisLoadAll(movie,dataProperties,idlist,loadOptions)
%IMARISLOADALL loads movies and idlists into imaris
%
% SYNOPSIS: imarisHandle = imarisLoadAll(movie,dataProperties,idlist,loadOptions)
%
% INPUT movie - 5D-array (x,y,z,c,t) or cell with one or several lines of
%               the form {fullMovieName,movieType,dimension}
%                   fullMovieName - fullpath of movie
%                   movieType - 'fim','stk','r3d'
%                   dimension - dimension along which to catenate movies.
%                               Optional if only 1 movie supplied.
%       dataProperties - dataProperties structure or cell of structures or
%                        filenames
%       idlist - optional idlist or cell of idlists. Group idlists in rows
%                that should be plotted into the same wavelength. The first
%                idlist will be plotted with the first wavelength, etc.
%       loadOptions - optional structure with fields
%                     .maxSize : maximum size of single data block
%                     .crop: cell with 5 elements (x,y,z,c,t) containing
%                         ranges of what to display
%
%
% OUTPUT imarisHandle - handle to the imaris application
%
% REMARKS
%
% created with MATLAB VERSION : 7.1.0.246 (R14) Service Pack 3 Windows_NT
%
% USERNAME: jdorn
% DATE: 27-Jan-2006
%
%

% defaults
maxSize = 100000000; % 100MB
crop = []; % no crop

%===================
% TEST INPUT
%===================

% check for required input
if nargin < 2 || isempty(movie) || isempty(dataProperties)
    error('please specify at least a movie with its corresponding dataProperties')
end

% check input individually
if iscell(movie)
    movieIsCell = 1;
    if size(movie,2) < 2
        error('movie input cell has to contain at least the filename and the movie type')
    end
    % remember number of movies
    nMovies = size(movie,1);
    % store movie information as movieCell
    movieCell = movie;
else
    movieIsCell = 0;
    movieCell = [];
end

% dataProperties. Make cell if movie is a cell. Don't check whether it is a
% struct or not.
if iscell(dataProperties)
    if length(dataProperties) ~= nMovies
        if length(dataProperties) == 1
            dataPropertiesCell = cell(nMovies,1);
            dataPropertiesCell{:} = deal(dataProperties{1});
        else
            error('number of movies and number of dataProperties is not equal')
        end
    end
else
    if movieIsCell
        dataPropertiesCell = cell(nMovies,1);
        dataPropertiesCell{:} = deal(dataProperties);
    end
end

% idlist. Only check if there is something. We check number of colors etc.
% later.
if isempty(idlist)
    isIdlist = 0;
else
    isIdlist = 1;
end

% loadOptions.
if ~isempty(loadOptions)
    if isfield(loadOptions,'maxSize')
        maxSize = loadOptions.maxSize;
    end
    if isfield(loadOptions,'crop')
        crop = loadOptions.crop;
    end
end

%==========================


%==========================
% PREPARE LOADING 
%==========================

% first of all, we need to calculate the new movieSize, and estimate which
% frames 
