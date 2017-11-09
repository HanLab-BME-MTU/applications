function [ sorted_fit ] = sortVonMisesFischer2d( fit )
%sortVonMisesFischer2d Sorts the peaks from a vonMisesFischer peak
%according to amplitude

if(iscell(fit))
    sorted_fit = cellfun(@sortVonMisesFischer2d,fit,'UniformOutput',false);
    return;
end

sorted_fit = fit;
% should be two columns
% first column: angle
% second column: amplitude
M = reshape(fit(1:end-2),2,[])';
% sort by amplitude then angle
M = flipud(sortrows(M(:,[ 2 1])));
sorted_fit(1:1:end-2) = reshape(M(:,[ 2 1])',1,[]);

end

