function [ xg_out, Kg ] = newtonBPprotoSimple( R, r, c)
%newtonBPprotoSimple Summary of this function goes here
%   Detailed explanation goes here

rho = R.getResponseAtOrderFTatPoint(r,c,R.filter.K);
xg = interpft_extrema(rho);
xg = xg(~isnan(xg));

% keyboard;
coords.r = repmat(r,1,length(xg));
coords.c = repmat(c,1,length(xg));
coords.m = 1:length(xg);
xg_out = zeros(1,length(xg));
Kg = zeros(1,length(xg));
for n = 1:length(xg)
    [xg_out(n), Kg(n)] = orientationSpace.diffusion.newtonBPproto(R, n, coords, xg(n), R.filter.K);
end

end

