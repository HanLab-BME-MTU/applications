function debugGroup

load(['C:\amatov\data\786O\786O_parental\786Opar_NaCl02_R3D\groups\group.mat']);

list = ([group.list]);
listU = unique(list);
leLu = length(listU);
k = 0;
for i = 1:leLu
   track = find(list==listU(i));
   if length(track)>1
       list(track)
       k = k + 1;      
       repTrack(k) = list(track(1));
   end
end
m = 0;
for j = 1:k
    for i = 1:length(group)
        aux = find(group(i).list==repTrack(j));
        if ~isempty(aux)
            m = m + 1;
            groupRepTrack(m,j) = i
        end
    end
end
group