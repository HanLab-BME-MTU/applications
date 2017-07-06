%   r = rand(7,1);
%   figure;
%   plot((0:length(r)-1)/length(r)*2*pi,r,'ko');
%   hold on;
%   plot((0:199)/200*2*pi,interpft(r,200),'k');
%   interpft_extrema(r);
%   hold off;

% r_hat = fft(r);
% half_length = floor(length(r)/2);
% f = [0:floor(half_length) -floor(half_length):-1].';
% 
% n_hat = r_hat;
% h = plot(r);
% ylim(ylim);
% for i=1:101
%     n = ifft([n_hat(1:half_length+1); zeros(354,1); n_hat(half_length+2:end)]);
%     h.YData = real(n(:).')*361/2/pi;
% %     plot(n);
%     n_hat = n_hat.*exp(-f.^2/2/half_length.^2);
%     pause(0.1);
%     hold on;
% end

%% Load image

if(~exist('I','var'))
    [~,hostname] = system('hostname');
    hostname = strtrim(hostname);

    switch(hostname)
        case 'FSM2UA220003Q'
            % HP Z820 Goldman workstation
            cd 'Z:\Takeshi\N-SIM\040715';
            MD = MovieData.load('MEFLB1-LACLB12-006_Reconstructed.nd2');
            I = MD.channels_(1).loadImage(1,11);
        case 'mkitti-jaqaman'
            % Laptop T440s
            load('C:\Users\Mark Kittisopikul\Documents\Data\Lamins\MEFLB1-LACLB12-006_Reconstructed_study\MEFLB1-LACLB12-006_Reconstructed\MEFLB1-LACLB12-006_Reconstructed.mat');
            MD.sanityCheck;
            I = MD.channels_(1).loadImage(1,10);
        case 'FSMPC0KTM9U'
            cd 'P:\Basic_Sciences\CMB\GoldmanLab\Takeshi\N-SIM\040715';
            MD = MovieData.load('MEFLB1-LACLB12-006_Reconstructed.nd2');
            I = MD.channels_(1).loadImage(1,11);
        otherwise
            % BioHPC
            cd ~/shortcuts/MEFLB1-LACLB12-006_Reconstructed/
            MD = MovieData.load('MEFLB1-LACLB12-006_Reconstructed.mat');
            I = MD.channels_(1).loadImage(1,10);
    end
end

%% Setup filter

% I = imread('example.tif');
F = OrientationSpaceFilter.constructByRadialOrder(1/2/pi./2,1,8,'none');
R = F*I;
K = 8:-0.1:1;
% t = linspace(1/(2*8+1).^2,1/(2+1).^2,1000);
% K = 1/2./sqrt(t)-1/2;
% rho = zeros(17,length(K));
% for i=1:length(K)
% %     out(:,i) = interpft_extrema(R.getResponseAtOrderFTatPoint(628,323,K(i)));
%       rho(:,i) = R.getResponseAtOrderFTatPoint(628,323,K(i));
% end
% rho = R.getResponseAtOrderFTatPoint(623,383,K);
% rho = R.getResponseAtOrderFTatPoint(628,323,K);
% rho = R.getResponseAtOrderFTatPoint(622,363,K);
r = 622;
c = 364;
rho = R.getResponseAtOrderFTatPoint(r,c,K);
out = interpft_extrema(rho);
out = orientationSpace.diffusion.alignExtrema(out);


outg = gradient(out)/0.1;
vderivs = interpft1_derivatives(rho,out,[1 2 3 4 5]);
% D = pi^2/2; % If we only using 180 degrees
D = 2*pi^2;
t = 1./(2*K+1).^2;
dt_dK = -4./(2*K+1).^3;
d2t_dK2 = 24./(2*K+1).^4;
dm_dt = -D*vderivs(:,:,3)./vderivs(:,:,2);
d2m_dt2 = -vderivs(:,:,3).*dm_dt.^2 - 2*D*vderivs(:,:,4).*dm_dt - D.^2.*vderivs(:,:,5);
d2m_dt2 = d2m_dt2 ./ vderivs(:,:,2);
% d2m_dK2 = d2m_dt2.*(dt_dK).^2 + dm_dt.*d2t_dK2;
d2m_dK2 = bsxfun(@times,d2m_dt2,dt_dK.^2) + bsxfun(@times,dm_dt,d2t_dK2);
% dm_dK = dt_dK.*dm_dt;
dm_dK = bsxfun(@times,dt_dK,dm_dt);


figure;
hold on;
plot(out.'/2/pi*180,interpft1([0 2*pi],rho,out).')
for i=0:7
    hold on; plot((0:359)/360*180,interpft(rho(:,1+10*i),360),'Color',ones(3,1)*0.125*i)
end
% xlim([0 2*pi]);
xlim([0 180]);
legend([strcat({'LM Track '},num2cell(num2str((1:size(out,1)).'))); strcat({'K = '},num2cell(num2str((8:-1:1).')))]);
grid on;

figure;
plot(K,out);
legend([strcat({'LM Track '},num2cell(num2str((1:size(out,1)).')))]);
grid on;

for trackNum = 1:size(out,1)
    try
        figure;
        title(['Track ' num2str(trackNum)]);
        idx_select = length(K):-35:1;
        track = out(trackNum,idx_select);
        idx_select = idx_select(~isnan(track));
        track = track(~isnan(track));

        x = joinColumns(repmat(K(idx_select),2,1));


        ygrad = joinColumns([track; -outg(trackNum,idx_select)]);
        spgrad = spapi(augknt(K(idx_select),4,2),x.',ygrad.');

        y = joinColumns([track; dm_dK(trackNum,idx_select)]);
        sp = spapi(augknt(K(idx_select),4,2),x.',y.');

        x2 = joinColumns(repmat(K(idx_select),3,1));
        y2 = joinColumns([track; dm_dK(trackNum,idx_select); d2m_dK2(trackNum,idx_select)]);
        sp2 = spapi(optknt(x2.',5),x2.',y2.');
        sp2c{trackNum} = sp2;


        xq = K(idx_select(1)):0.01:K(idx_select(end));
        xqc{trackNum} = xq;
                hold on;
        plot(xq,spval(sp,xq))
        plot(xq,spval(spgrad,xq))
        plot(xq,spval(sp2,xq),'--')
        plot(K,out(trackNum,:),'o');
        plot(K(idx_select),out(trackNum,idx_select),'o','MarkerFaceColor','k');
    catch err
    end
end

figure;
title('First Derivatives');
hold on;
plot(K,-outg)
hold on;
plot(K,dm_dK,'o')
grid on

figure;
title('Absolute Value First Derivatives');
% plot(K,abs(outg))
hold on;
plot(K,abs(dm_dK))
grid on

figure;
plot(K,d2m_dK2)
title('Second Derivative');
grid on
hold on;
plot(K,gradient(dm_dK,-0.1),'o')