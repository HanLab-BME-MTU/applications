function [ varargout ] = demo_nlms_2( complex, radius , noise)
%demo_nlms_2 Demonstrate nlms

if(nargin < 1 || isempty(complex))
    complex = true;
end
if(nargin < 2 || isempty(radius))
    radius = 35;
end
if(nargin < 3 || isempty(noise))
    noise = 1e-3;
end

% I = zeros(101,101);
[X,Y] = meshgrid(-50:50,-50:50);
R = hypot(X,Y);
I = R > radius-1/2 & R < radius+1/2;
if(complex)
    I = I | circshift(I,51,1) | circshift(I,51,2);
end
I = double(I);
I = imfilter(I,fspecial('gaussian',10,2));
I = imnoise(I,'gaussian',0,noise);

[v.res, v.theta, v.nms, v.a] = steerableOrientationSpaceFilter(I,0.08,[],5,180);
% [v.res, v.theta, v.nms, v.a] = steerableDetector(I,4,5,256);
nlms = nonLocalMaximaSuppression(real(v.a));
T = thresholdOtsu(nlms(nlms > 0));
spy3d(nlms > T)
hlines = findobj(gca,'Type','Line');
hold on;
hSurf = surf(0.5:101,0.5:101,-ones(101),imfuse(I,real(v.nms)),'EdgeColor','none');
zlim([-1 180]);
colormap gray
if(nargout > 0)
    varargout{1} = nlms;
end

hSpin = uicontrol('Style','pushbutton','String','Spin','Callback',@spin);
hToggle = uicontrol('Style','togglebutton','String','Toggle Spy','Callback',@toggleSpy);
pos = get(hToggle,'Position');
set(hToggle,'Position',pos+[100 0 0 0]);

    function spin(source,cbdata)
        phiDir = 1;
        for i=1:360;
            [az,el] = view;
            camorbit(1,phiDir);
            if(el < 0)
                phiDir = -phiDir;
                camorbit(0,phiDir*2);
            end
            disp(el);
            disp(phiDir);
            drawnow;
%             pause(0.05);
        end
    end
    function toggleSpy(source,cbdata)
        visState = get(hlines(1),'Visible');
        switch(visState)
            case 'on'
                set(hlines,'Visible','off');
            case 'off'
                set(hlines,'Visible','on');
        end
    end


end

