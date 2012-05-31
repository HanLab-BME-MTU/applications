function plotLifetimesEAF(data)

t = 1:data.movieLength;

figure;
ni = hist([data.tracks.lifetime], t);
ni = ni / sum(ni);
plot(t, ni, 'k.-');
hold on;

idx = [data.tracks.cStatus];
idx = idx(2,:);

ni = hist([data.tracks(idx==1).lifetime], t);
ni = ni / sum(ni);
plot(t, ni, 'g.-');

ni = hist([data.tracks(idx==0).lifetime], t);
ni = ni / sum(ni);
plot(t, ni, 'r.-');

legend('All CCPs', 'EAF+', 'EAF-');