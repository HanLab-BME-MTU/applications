function job = makiMakeJobPlatformIndependent(job,serverType)
%MAKIMAKEJOBPLATFORMINDEPENDENT makes a job platform independent
%
% SYNOPSIS: job = makiMakeJobPlatformIndependent(job,serverType)
%
% INPUT job: maki job file
%       serverType: 'TEST','HERCULES','DANUSER','MERALDI',SWEDLOW' or
%       'MCAINSH'
%
% OUTPUT 
%
% REMARKS
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn
% DATE: 05-Jul-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% list here all fields that contain path names in dataStruct and job
dataFields2check = {'rawMoviePath';'dataFilePath'};
jobFields2check = {'jobPath'};

nDataFields = length(dataFields2check);
nJobFields = length(jobFields2check);

nJobs = length(job);

% go through all the fields of a job file and switch paths to identifiers
% or back
for iJob = 1:nJobs
   % check job fields
   for f=1:nJobFields
       job(iJob).(jobFields2check{f}) = ...
           makiPathDef(job(iJob).(jobFields2check{f}),serverType);
   end
   % check dataStruct fields
   for f=1:nDataFields
       job(iJob).dataStruct.(dataFields2check{f}) = ...
           makiPathDef(job(iJob).dataStruct.(dataFields2check{f}),serverType);
   end
    
end