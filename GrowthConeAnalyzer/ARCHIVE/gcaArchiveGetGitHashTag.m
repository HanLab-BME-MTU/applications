function [success,hashTag] = gcaArchiveGetGitHashTag
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

OS = getenv('OS') ;
if ~isempty(regexp(OS,'Windows','once'));
    load('C:\Users\Maria\matlab\REPOSITORY_GIT\applications\GrowthConeAnalyzer\hashTag.mat');
    % i am running it - later figure out how to get the hash tag
    % generically...
else
    % try to load hashtag (likely on linux)
    cDir = pwd;
    cd(['/home2' filesep 'mbagonis' filesep 'matlab' filesep 'applications']);
    cmd = ['TERM=ansi git log --pretty=%H -2'];
    [success,hashTag] =  system(cmd);
    cd(cDir);
    save([cDir filesep 'hashTag.mat'],'hashTag'); 
end



end

