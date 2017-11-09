function [ bifurcation_points ] = findLocalMaxExtremum( response, K )
%findBifurcation Find bifurcation points

response_d = ifft(fft(response).*ifftshift((-K(1):K(1))*1i).');
response_dd = ifft(fft(response).*ifftshift((-K(1):K(1))*1i).'.^2);
% Just to get the expanded repsonse_d
[spsmax,spsmin,response_d] = orientationSpace.diffusion.splinesForExtremaFT(response_d,K);
[spsmax,spsmin,response_dd] = orientationSpace.diffusion.splinesForExtremaFT(response_dd,K);
sps = [spsmax spsmin];

% validK = arrayfun(@(x) K <= max(x.knots) & K >= min(x.knots),[spsmax spsmin],'Unif',false);
% test = arrayfun(@(x,validK) interpft1([0 2*pi],rhod(:,validK{1}),spval(x,K(validK{1}))),[spsmax spsmin],validK,'Unif',false)

firstDerivZeros  = cell(1,length(sps));
firstDerivZeroValues = cell(1,length(sps));
for s = 1:length(sps)
    validK = K <= max(sps(s).knots) & K >= min(sps(s).knots);
    firstDerivValues = interpft1([0 2*pi],response_d(:,validK),spval(sps(s),K(validK)));
    firstDerivSpline = spline(K(validK),firstDerivValues);
    firstDerivZeros{s} = fnzeros(firstDerivSpline);
    firstDerivZeros{s} = firstDerivZeros{s}(1,firstDerivZeros{s}(1,:) == firstDerivZeros{s}(2,:));
    if(~isempty(firstDerivZeros{s}))
        firstDerivZeroValues{s} = spval(sps(s),firstDerivZeros{s}(1,:));
    end
end

bifurcation_points = [firstDerivZeros{:}; firstDerivZeroValues{:}].';

end

