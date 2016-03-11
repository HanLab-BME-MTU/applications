function [ h ] = vanGinkelVisualizeDeletion( orientationMatrix, target, r, c)
%vanGinkelUpsample upsample orientation information
%
% INPUT
% orientationMatrix : YxXxOrientation
% target            : Orientation angles to estimate
%
% OUTPUT
% y                 : YxXxTarget
    import vanGinkel.*;
    if(nargin < 2)
        target = pi/72;
    end
    
    % 
    s = size(orientationMatrix);
    
    if(nargin < 3)
        r = round(s(1)/2);
    end
    
    if(nargin < 4)
        c = round(s(2)/2);
    end
    
    ind = sub2ind(s(1:2),r,c);
    
    M = reshape(orientationMatrix,s(1)*s(2),s(3));
    % todo, deal with even case
    n = (s(3))/2;
%     n = n+0.5;

    
%     x = [0:n -n:-1];
    x = wraparoundN(0:s(3)-1,-n,n);
    % toeplitz?
    xx = wraparoundN(bsxfun(@minus,x,(0:(2*n-1))'),[-n n]);
    % xx = arrayfun(@(k) circshift(x,k,2),0:s(3),'UniformOutput',false);
    % xx = vertcat(xx{:});

    A = 2*exp(-xx.^2/2);
    
    [vderiv,vderiv2] = vanGinkelDerivative(orientationMatrix,target);


    if(isscalar(target))
        target = 0:target:pi-target;
    end
    target = target./(pi/s(3));
    
    % note we factor out the exp(1/2) constant factor

    T = bsxfun(@minus,x,target')';
    T = wraparoundN(T,[-n n]);
    T = 2*exp(-T.^2/2);
    w = A\T;
    % normalize the weights such that the columns sum to 1
%     w = w./repmat(sum(w),s(3),1);

    y = M*w;
%     y = reshape(y,s(1),s(2),[]);

    
    % Visualize
    plots = @(f,varargin) plot((0:length(f)-1)/length(f),squeeze(real(f)),varargin{:});
    
    h.fig = figure;
%     h.mainCircles = plots(M(ind,:),'ko');
    hold on;
    h.mainLine = plots(y(ind,:),'k','LineWidth',2);
    h.gaussians = plots(bsxfun(@times,(M(ind,:)/A)',T));
    set(h.gaussians,'ButtonDownFcn',@toggleLine);
%     h.bessels = plots(bsxfun(@times,M(ind,:)',A\T));
%     set(h.bessels,'Visible','off');
%     plots(bsxfun(@times,M(ind,:)',w));
%     h.derivatives(1) = plots(vderiv(r,c,:),'k--');
%     h.derivatives(2) = plots(vderiv2(r,c,:),'k:');
    h.ax = gca;
    grid on;
    set(gca,'XTick',(0:s(3)-1)/s(3))
    set(gca,'XTickLabel',strcat(cellfun(@num2str,num2cell(0:s(3)-1),'Unif',false),'/',num2str(s(3))));
    
%     pos = get(h.ax,'Position');
%     pos(4) = 0.8 - pos(2);
    set(h.ax,'Position',[0.05 0.1 0.9 0.85],'XLimMode','manual','YLimMode','manual','XLim',xlim,'YLim',ylim);
% %     pos(2) = 0.8;
%     h.panel = uipanel(h.fig,'Position',[0.05 0 0.9 0.05]);
%     pos = [0 0 100 20];
%     h.toggles.derivatives = uicontrol(h.panel,'Style','togglebutton','String','Derivatives','Position',pos,'Value',1, ...
%         'BackgroundColor',[0 1 0],'Callback',@toggleCallback,'UserData','derivatives');
%     pos(1) = pos(1) + 120;
%     h.toggles.gaussians = uicontrol(h.panel,'Style','togglebutton','String','Gaussians','Position',pos,'Value',1, ...
%         'BackgroundColor',[0 1 0],'Callback',@toggleCallback,'UserData','gaussians');
%     pos(1) = pos(1) + 120;
%     h.toggles.bessels = uicontrol(h.panel,'Style','togglebutton','String','Bessels','Position',pos,'Value',0, ...
%         'BackgroundColor',[1 0 0],'Callback',@toggleCallback,'UserData','bessels');
%     pos(1) = pos(1) + 120;
%     h.toggles.oneBessel = uicontrol(h.panel,'Style','pushbutton','String','One Bessel','Position',pos,'Value',0, ...
%         'Callback',@oneBesselVisible,'UserData','bessels');
    guidata(h.fig,h);
    

end
function toggleLine(hline,ev)
    h = guidata(hline);
    grey = [0.9 0.9 0.9];
    if(all(hline.Color == grey))
        % restore gaussian
        hline.Color = hline.UserData;
        h.mainLine.YData = h.mainLine.YData + hline.YData;
    else
        % delete gaussian
        hline.UserData = hline.Color;
        hline.Color = grey;
        h.mainLine.YData = h.mainLine.YData - hline.YData;
    end
end
% function toggleCallback(hButton,ev)
%     h = guidata(hButton);
%     states = {'off','on'};
%     colors = {[1 0 0],[0 1 0]};
%     set(h.(hButton.UserData),'Visible',states{hButton.Value+1});
% %     switch(hButton.String)
% %         case 'Derivatives'
% %             set(h.(hButton.UserData),'Visible',states{hButton.Value+1});
% %         case 'Gaussians'
% %             set(h.(hButton.UserData),'Visible',states{hButton.Value+1});
% %     end
%     set(hButton,'BackgroundColor',colors{hButton.Value+1});
% end
% function oneBesselVisible(hButton,ev)
%     h = guidata(hButton);
%     set(h.bessels,'Visible','off');
%     set(h.bessels(randi(length(h.bessels))),'Visible','on');
% end