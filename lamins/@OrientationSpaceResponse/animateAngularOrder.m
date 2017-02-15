function [  ] = animateAngularOrder( R, r, c, Rd )
%animateAngularOrder Summary of this function goes here
%   Detailed explanation goes here

if(nargin < 2)
    r = 127;
end
if(nargin < 3)
    c = 126;
end
if(nargin < 4)
    da = ifft(bsxfun(@times,fft(R.a,[],3),shiftdim([0:8 -8:-1],-1))*1i,[],3);
    Rd = OrientationSpaceResponse(R.filter,da);
end

pauseTime = 0.1;
deltaK = 0.1;

h(1) = plot(0:0.5:179.5,interpft(R.getResponseAtOrderFTatPoint(r,c,8),360));
hold on;
h(2) = plot(0:0.5:179.5,interpft(Rd.getResponseAtOrderFTatPoint(r,c,8),360));
h(3) = line;
h(3).LineStyle = 'none';
h(3).Marker = '.';
h(3).MarkerSize = 1;
h(4) = line;
h(4).Marker = 'o';
h(4).LineStyle = 'none';
grid on;
drawnow;
yl = ylim;
xl = xlim;
ht = text(xl(2)-30,yl(2)-10,['K = ' num2str(R.filter.K)]);
ylim('manual');


while all(isvalid(h))
    h(3).XData = [];
    h(3).YData = [];
    for K=R.filter.K:-deltaK:-0.5+deltaK
        temp = R.getResponseAtOrderFTatPoint(r,c,K);
        h(1).YData = interpft(temp,360);
        h(2).YData = interpft(Rd.getResponseAtOrderFTatPoint(r,c,K),360);
        ht.String = ['K = ' num2str(K)];
        [maxima,~,maxima_value] = interpft_extrema(R.getResponseAtOrderFTatPoint(r,c,K));
        h(3).XData = [h(3).XData maxima.'/2/pi*180];
        h(3).YData = [h(3).YData maxima_value.'];
        h(4).XData = maxima.'/2/pi*180;
        h(4).YData = maxima_value.';
        pause(pauseTime);
    end
end

end

