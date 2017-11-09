function [ h ] = scatter_prctile( X, Y , pvec , varargin)
%scatter_prctile Scatters of the percentiles indicated by pvec of samples X
%and Y

    Xp = prctile(X,pvec);
    Yp = prctile(Y,pvec);
    h = scatter(Xp,Yp,varargin{:});
    allmax = max([Xp Yp]);
    ylim([0 allmax]);
    xlim([0 allmax]);

end

