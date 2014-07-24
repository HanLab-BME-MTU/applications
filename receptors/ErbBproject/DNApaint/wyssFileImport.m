function [PointList] = wyssFileMeanShift(dir,varargin)
% takes in a directory that contains .mat files with data from the wyss
% converts the data into a PointList structure and then runs meanshift
% clustering on the data.
%
% fullldata has the follow columns
% (0) frame
% (1) x
% (2) y
% (3) photons
% (4) sigma_x
% (5) sigma_y
% (6) loc_ac
% (7) offset
% (8) noise
% (9) amp
%
% Results are then saved in the parent directory of dir as "dir".mat
%

ip = inputParser;

ip.addRequired('dir',@ischar);

ip.addOptional('name',[],@(x) ischar(x)|isempty(x));

ip.parse(dir,varargin{:});

name = ip.Results.name;

cd(dir);
if strcmp(dir(end),filesep)
    dir = dir(1:end-1);
end
output = [dir,'.mat'];

list = what;
list = list.mat;


PointList = cell([numel(list),1]);
TotalPnts = [];

% determines which paint imaging strand was used and reorder the point list
% so that it starts with p1 and end with p13. assumes files name starts
% with pXX
ind = [4,5,1,2,3];

for i = 1:numel(list);
    
load(list{i});
PointList(ind(i)) = {struct('pnts',trace(:,2:3),'name',list{i},'fullData',trace,'drift',[],'dmark',[],'shift',[],'com',[])};
TotalPnts = vertcat(TotalPnts,[trace(:,2:3),ind(i)*ones(size(trace(:,1)))]);

end

if ~isempty(name)
save([name,'_PointList.m'],'PointList');
end

end
