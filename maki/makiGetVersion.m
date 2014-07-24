function [versionNum,fileBody] = makiGetVersion(fName)
%MAKIGETVERSION fetches the version number of a file
%
%SYNOPSIS [versionNum,fileBody] = makiGetVersion(fName)
%
%INPUT  fName     : File name.
%
%OUTPUT versionNum: Version number.
%       fileBody  : File name body, without the underscore & version number. 
%
%The function assumes that the version index is the last number before
%'.mat' and is preceeded by an '_';
%
%for example: 'gaga_cpi11_anythingelse_23.mat'
%will generate a numerical value 23

indx = regexp(fName,'_\d+\.mat');

versionNum = str2num(fName(indx+1:regexp(fName,'.mat')-1));

fileBody = fName(1:indx-1);


            
