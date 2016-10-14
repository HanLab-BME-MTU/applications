function [ x,y,z] = upDirectory( dir,num,fp )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin<2 || isempty(num)
   num =1;  
end

if nargin<3 || isempty(fp)
    fp = 0 ; 
end 
x = dir; 


for i = 1:num
    if fp == 0
   [x,y,z] =  getFilenameBody(x); 
    else 
        [x,y,z] = fileparts(x); 
    end 
end

    