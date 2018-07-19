function [ x,y,z] = mitoticUpDirectory( dir,num,fp )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin<3 || isempty(fp)
    fp = 0 ; 
end 
x = dir; 


for i = 1:num
    if fp == 0
   [x,y,z] =  mitoticGetFilenameBody(x); 
    else 
        [x,y,z] = fileparts(x); 
    end 
end

    