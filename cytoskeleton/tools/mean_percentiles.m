function MP = mean_percentiles(this_data)
this_data = this_data(:);
MP = nan(1,8);
% if there is nothing to get percentile, return nan
if( isempty(this_data) || sum(isnan(this_data))== numel(this_data))
    return;
end
MP(1) = nanmean(this_data(:));
MP(2) = prctile(this_data(:),0);
MP(3) = prctile(this_data(:),2);
MP(4) = prctile(this_data(:),25);
MP(5) = prctile(this_data(:),50);
MP(6) = prctile(this_data(:),75);
MP(7) = prctile(this_data(:),98);
MP(8) = prctile(this_data(:),100);
