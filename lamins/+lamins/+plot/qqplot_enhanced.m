function [ h ] = qqplot_enhanced( varargin )
%qqplot_enhanced Summary of this function goes here
%   Detailed explanation goes here

if(nargin > 2)
    quantiles = varargin{3};
else
    quantiles = 100 - [(sqrt(2).^-(0:50))*100 0];
end

A = varargin{1};
B = varargin{2};

% axes;
% hold on;

if(nargin > 3)
    A_labels = prctile(A,varargin{4});
    B_labels = prctile(B,varargin{4});
    arrayfun(@(A,B) line([A A 0],[0 B B],'Color',[0.75 0.75 0.75],'LineStyle','-'),A_labels,B_labels);
%     labels = round(varargin{4},1);
%     labels = arrayfun(@num2str,labels,'Unif',false);
%     text(A_labels,B_labels,labels);
end

A_med = nanmedian(A(:));
B_med = nanmedian(B(:));

h.median = line([A_med A_med 0],[0 B_med B_med],'Color',[0 0 0]);

hold on;

h.original = qqplot(varargin{1:2},quantiles);





axis square;

current.ylim = ylim;
ylim([0 current.ylim(2)]);



h.main_diag = line([0 current.ylim(2)],[0 current.ylim(2)],'Color',[0.5 0.5 0.5],'LineStyle',':');

h.lim = [0 current.ylim(2)];

h.scaleFactor = lsqr(prctile(A,25:75)',prctile(B,25:75)');

delete(h.original(2:3));

h.original(2) = line([0 h.lim(2)],[0 h.lim(2).*h.scaleFactor],'Color','r','LineStyle','-.');

h.ranksum = ranksum(A*h.scaleFactor,B);

hold off;



end

