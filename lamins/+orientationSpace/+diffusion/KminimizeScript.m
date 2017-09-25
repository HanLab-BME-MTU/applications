maxima = R.getRidgeOrientationLocalMaxima;
import orientationSpace.diffusion.*;
[ maxima, coords, response] = linearizeMaxima( R, maxima);
[ d.maxima, d.coords, d.response, d.xgd ] = descendHalleyK( maxima, coords, response );
[x_out,K_out,notDone] = refineBifurcation(d.maxima,d.coords.K,response,8);
maximaF = NaN(1024,1024,7);
coords.ind = sub2ind([1024 1024 7],coords.r,coords.c,coords.m);
maximaF(coords.ind) = x_out;
[fr,fc] = find(sum(maximaF(:,:,1:end-1) == maximaF(:,:,2:end),3));
maximaF_repeats = cat(3,false(1024),maximaF(:,:,1:end-1) == maximaF(:,:,2:end));
nn = find(maximaF_repeats(coords.ind));
x_fixed = x_out;
K_fixed = K_out;
for n = nn
[x_fixed(n),K_fixed(n)] = orientationSpace.diffusion.newtonBPproto(R, n, coords, maxima(n), 8);
end
maximaFixed = NaN(1024,1024,7);
maximaFixed(coords.ind) = x_fixed;
maximaF_repeats = cat(3,maximaF(:,:,1:end-1) == maximaF(:,:,2:end),false(1024));
nn = find(maximaF_repeats(coords.ind));
for n = nn
[x_fixed(n),K_fixed(n)] = orientationSpace.diffusion.newtonBPproto(R, n, coords, maxima(n), 8);
end
save('maximaFixed.mat','maximaFixed')

K_fixedN = NaN(1024,1024,7);
K_fixedN(coords.ind) = K_fixed;
maximaFixedN(K_fixedN <= 3) = NaN;
maximaCombined = cat(3,maximaFixedN,maxima3);
maximaCombined = sort(maximaCombined,3);
max(max(sum(~isnan(maximaCombined),3)))
maximaCombined = maximaCombined(:,:,1:7);