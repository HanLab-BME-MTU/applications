function [ output_args ] = GCAVisualsProtrusionMapMovie(movieData)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% INPUT
%
%   movieData - The MovieData object describing the movie, as created using
%   movieSelectorGUI.m
%
%   Parameter Structure Field Names:
%
% Generic Fields: (Input/Output Fields Needed for Wrapper)
%       ('OutputDirectory' -> Optional. A character
%       string specifying the directory to save the Visualization Output
% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end

% load the protrusionSamples 

%% Check Input 
cDir = [MD.outputDirectory_ filesep 'protrusion_samples'];
load([cDir filesep 'protrusion_samples.mat']); 

GCAVisualsProtrusionMap(protSamples,cDir,0); 
end

