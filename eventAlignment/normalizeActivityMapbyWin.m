function normActivity = normalizeActivityMapbyWin(activity)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[nWin nTime]=size(activity);
for j=1:nWin
    minV=min(activity(j,:));
    maxV=max(activity(j,:));
    normActivity(j,:)=(activity(j,:)-minV)/(maxV-minV)*1000;
end


