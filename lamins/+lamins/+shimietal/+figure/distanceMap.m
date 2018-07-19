function [ data, dataMedian ] = distanceMap( stats )
%distanceMap
%
% merged structture of lamin data
% e.g. load('/project/biophysics/jaqaman_lab/lamins/2015/20150602/stats_20160624.mat')

% labels = {'LAC', 'LB1', 'LB2'};
% 
% WT.idx = 20;
% % channels for comparison
% % Use LAC channel for LB nulls
% WT.LB1 = 1;
% WT.LB2 = 1;
% % Use LB12 channel for LAC nulls
% WT.LAC = 2;
% 
% idx = [2 7 12];
% 
% immuno.A.idx = 21;
% immuno.B1.idx = 21;
% immuno.B2.idx = 28;
% immuno.C.idx = 28;
% 
% immuno.A.ch = 1;
% immuno.B1.ch = 2;
% immuno.B2.ch = 2;
% immuno.C.ch = 1;
% 
% mEmerald.A.idx = 29;
% mEmerald.B1.idx = 30;
% mEmerald.B2.idx = 31;
% mEmerald.C.idx = 32;
% 
% mEmerald.A.ch = 1;
% mEmerald.B1.ch = 1;
% mEmerald.B2.ch = 1;
% mEmerald.C.ch = 1;
% 
% idx = [21 22 25 27 28 23]';
%                 
% hist.labels = {  ...
%    'LA'  'LB1';  ...
%    'LA'  'LB2';  ...
%    'LB2' 'LB1';  ...
%    'LC'  'LB1';  ...
%    'LC'  'LB2';  ...
%    'LC'  'LA' };

% label idx ch
masterData = {   ...
    '\alpha-LA'    21 1 ;...
    '\alpha-LB1'   21 2 ;...
    '\alpha-LB2'   22 2 ;...
    '\alpha-LC'    23 1 ;...
    'mLA'     29 1 ;...
    'mLB1'    30 1 ;...
    'mLB2'    31 1 ;...
    'mLC'     32 1 ;...
    '\alpha-LAC'   20 1 ;...
    '\alpha-LB12'  20 2 ;...
    '\alpha-LB12 LA-/-'   2 1 ;...
    '\alpha-LAC LB1-/-'  7 1 ;...  
    '\alpha-LAC LB2-/-' 12 1
    };

%  1. a-LA
%  2. a-LB1
%  3. a-LB2
%  4. a-LC
%  5. mLA
%  6. mLB1
%  7. mLB2
%  8. mLC
%  9. a-LAC
% 10. a-LB12
% 11. a-LB12 LACnull
% 12. a-LAC LB1null
% 13. a-LAC LB2null

c = [1 1 1 1 2 2 2 2 3 3 3 3 3];
% c = 'rrrrggggbbbbb';

fields = {'edgesPerFace','avgEdgeLength'};

for f=1:2
    field = fields{f};
    for i=1:size(masterData,1)
        for j=1:size(masterData,1)
            Xi.idx = masterData{i,2};
            Xi.ch = masterData{i,3};
            X = stats(Xi.idx).(field)(Xi.ch).all.data;
            Yi.idx = masterData{j,2};
            Yi.ch = masterData{j,3};
            Y = stats(Yi.idx).(field)(Yi.ch).all.data;
            data.(field)(i,j) = lamins.functions.qqscale(X,Y);
        end
    end
%     distM.(field) = ( data.(field) + data.(field)')/2;
end

for f=1:2
    field = fields{f};
    for i=1:size(masterData,1)
            Xi.idx = masterData{i,2};
            Xi.ch = masterData{i,3};
            X = nanmedian(stats(Xi.idx).(field)(Xi.ch).all.data);
            dataMedian.(field)(i) = X;
    end
%     distM.(field) = ( data.(field) + data.(field)')/2;
end

Twt = wtTable(stats);

data.edgesPerFace(1,3) = Twt.edgesPerFace(2);
data.edgesPerFace(1,4) = 1./Twt.edgesPerFace(6);

data.avgEdgeLength(1,3) = Twt.avgEdgeLength(2);
data.avgEdgeLength(1,4) = 1./Twt.avgEdgeLength(6);

% scaling these changes by the LA median
xScale = dataMedian(1,1).edgesPerFace(1,1);
yScale = dataMedian(1,1).avgEdgeLength(1,1);

figure;
scatter(data.edgesPerFace(1,c == 1).*xScale,data.avgEdgeLength(1,c == 1).*yScale,[],'r','o');
hold on;
scatter(data.edgesPerFace(1,c == 2).*xScale,data.avgEdgeLength(1,c == 2).*yScale,[],'g','s');
scatter(data.edgesPerFace(1,c == 3).*xScale,data.avgEdgeLength(1,c == 3).*yScale,[],'b','x');
line([0.96 1.1].*xScale,[1 1].*yScale,'LineStyle','--','Color',[0.5 0.5 0.5]);
line([1 1].*xScale,[0.96 1.1].*yScale,'LineStyle','--','Color',[0.5 0.5 0.5]);
xlim([0.96 1.1].*xScale); ylim([0.96 1.1].*yScale);
text((data.edgesPerFace(1,:)+0.002).*xScale,(data.avgEdgeLength(1,:)).*yScale,masterData(:,1));
xlabel('Number of Edges per Face (scaled)');
ylabel('Average Edge Length per Face (scaled, \mum)');
% xlabel('Edges per Face, fold change over \alpha-LA');
% ylabel('Average Edge Length, fold change over \alpha-LA');

% line_color = [0.75 0.75 0.75];
% line(data.edgesPerFace(1,[1 2]),data.avgEdgeLength(1,[1 2]),'Color',line_color);
% line(data.edgesPerFace(1,[1 3]),data.avgEdgeLength(1,[1 3]),'Color',line_color);
% line(data.edgesPerFace(1,[1 4]),data.avgEdgeLength(1,[1 4]),'Color',line_color);
% 
% line(data.edgesPerFace(1,[1 5]),data.avgEdgeLength(1,[1 5]),'LineStyle',':','Color',line_color);
% line(data.edgesPerFace(1,[2 6]),data.avgEdgeLength(1,[2 6]),'LineStyle',':','Color',line_color);
% line(data.edgesPerFace(1,[3 7]),data.avgEdgeLength(1,[3 7]),'LineStyle',':','Color',line_color);
% line(data.edgesPerFace(1,[4 8]),data.avgEdgeLength(1,[4 8]),'LineStyle',':','Color',line_color);
% 
% line(data.edgesPerFace(1,[9 12]),data.avgEdgeLength(1,[9 12]),'LineStyle','--','Color',line_color);
% line(data.edgesPerFace(1,[9 13]),data.avgEdgeLength(1,[9 13]),'LineStyle','--','Color',line_color);
% 
% line(data.edgesPerFace(1,[10 11]),data.avgEdgeLength(1,[10 11]),'LineStyle','--','Color',line_color);

colormap(eye(3));




% 
% figure;
% scatter(log(data.edgesPerFace(1,:)),log(data.avgEdgeLength(1,:)),[],c);
% line([-0.04 0.1],[0 0],'LineStyle','--','Color',[0.7 0.7 0.7]);
% line([0 0],[-0.04 0.1],'LineStyle','--','Color',[0.7 0.7 0.7]);
% xlim([-0.04 0.1]); ylim([-0.04 0.1]);
% text(log(data.edgesPerFace(1,:)+0.002),log(data.avgEdgeLength(1,:)),masterData(:,1));
% xlabel('Edges per Face, log ratio over \alpha-LA');
% ylabel('Average Edge Length, log ratio over \alpha-LA');
% colormap(eye(3));




   

end

