function [ bifurcation_points ] = findBifurcation( response, K )
%findBifurcation Find bifurcation points

response_d = ifft(fft(response).*ifftshift((-K(1):K(1))*1i).');
[spsmax,spsmin,response_d] = orientationSpace.diffusion.splinesForExtremaFT(response_d,K);
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

