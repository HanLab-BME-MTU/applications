function imagesc_white_nan(data, minValue, MaxValue)

if(nargin<2)
    minValue = min(min(data(~isnan(data))));
end
if(nargin<3)
    MaxValue = max(max(data(~isnan(data))));
end

 [nr,nc] = size(data);
 pcolor([data nan(nr,1); nan(1,nc+1)]); shading flat; axis off;
  caxis([minValue MaxValue]);
