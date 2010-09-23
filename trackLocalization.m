

% new function that performs Gaussian fitting, saves to Detection/localizationResults.mat
function trackLocalization(data, sigma)


ny = data.imagesize(1);
nx = data.imagesize(2);

frameList = dir([data.source '*.tif*']);
maskPath = [data.source 'Detection' filesep 'Masks' filesep];
maskList = dir([maskPath '*.tif']);

w = ceil(3*sigma);


load([data.source 'Detection' filesep 'detectionResults.mat']);

% loop through frames
fprintf('Progress:     ');
for k = 1:data.movieLength
    
    frame = double(imread([data.source frameList(k).name]));
    mask = double(imread([maskPath maskList(k).name]));
    
    xi = round(frameInfo(k).xcom);
    yi = round(frameInfo(k).ycom);
    
    np = length(xi);
    xloc = zeros(1,np);
    yloc = zeros(1,np);
    I = zeros(1,np);
    cVect = zeros(1,np);
    
    for p = 1:np
        
        %xi = round(frameInfo(k).xcom(p));
        %yi = round(frameInfo(k).ycom(p));
        
        % detections within frame bounds
        if xi(p)<=w || xi(p)>nx-w || yi(p)<=w || yi(p)>ny-w
            xloc(p) = NaN;
            yloc(p) = NaN;
            I(p) = NaN;
            c(p) = NaN;
        else
            
            window = frame(yi(p)-w:yi(p)+w, xi(p)-w:xi(p)+w);
            %figure; imagesc(window); colormap(gray(256)); axis image;
            
            
            % binary mask
            maskWindow = ~mask(yi(p)-w:yi(p)+w, xi(p)-w:xi(p)+w);
            maskWindow(maskWindow~=0) = 1;
            % background estimate
            c = mean(mean(window(maskWindow)));
            prm = fitGaussian2D(window, [0 0 max(window(:))-c sigma c], 'xyA');           
            xloc(p) = xi(p) + prm(1);
            yloc(p) = yi(p) + prm(2);           
            I(p) = prm(3)*2*pi*sigma^2;
            cVect(p) = c;
        end
    end
    frameInfo(k).xloc = xloc;
    frameInfo(k).yloc = yloc;
    frameInfo(k).I = I;
    frameInfo(k).c = cVect;
    fprintf('\b\b\b\b%3d%%', round(100*k/(data.movieLength)));
end
fprintf('\n');

save([data.source 'Detection' filesep 'detectionResults.mat'], 'frameInfo');

