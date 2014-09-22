%

for iFrame = 1 : MD.nFrames_;
    
    I = double(MD.channels_.loadImage(iFrame));
    max_I = max(max(I));
    I = I/max(max(I));
    
    [XI, YI] = meshgrid(1:size(I,2), 1:size(I,1));
    
    fit_sur = fit([YI(:),XI(:)],I(:), 'poly22', 'Robust', 'on');
    
    Z_fit = fit_sur.p01.*XI+ fit_sur.p10.*YI +fit_sur.p11.*XI.*YI+ ...
        fit_sur.p20.*YI.*YI+fit_sur.p02.*XI.*XI;
    
    figure(1);imshow([I-Z_fit]);
    imwrite(uint16([I-Z_fit]*max_I), ['background_removed_',num2str(iFrame),'.tif']);
end