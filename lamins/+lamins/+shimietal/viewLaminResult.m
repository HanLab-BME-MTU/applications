% import lamins.functions.*;
% 
% %% Load Lamin A label
% LA_label.WT.ML = MovieList.load('MEFWT-LA/MEFWT-LA_list.mat');
% LA_label.LB1null.ML = MovieList.load('MEFLB1-LA/MEFLB1-LA_list.mat');
% LA_label.LB2null.ML = MovieList.load('MEFLB2-LA/MEFLB2-LA_list.mat');
% LA_label.mEmerald.ML = MovieList.load('MEFWTmEmerald-LAfixed/MEFWTmEmerald-LAfixed_list.mat');
% 
% %% Aggregate data
% LA_label.WT.data = lamins.functions.aggregateMovieListData(LA_label.WT.ML);
% LA_label.LB1null.data = lamins.functions.aggregateMovieListData(LA_label.LB1null.ML);
% LA_label.LB2null.data = lamins.functions.aggregateMovieListData(LA_label.LB2null.ML);
% LA_label.mEmerald.data = lamins.functions.aggregateMovieListData(LA_label.mEmerald.ML);
% 
% %% Get Statistics
% % [LA_label.WT.h, LA_label.WT.rps, LA_label.WT.all, LA_label.WT.areas] = plotFaceStatistics([LA_label.WT.data.Sf]);
% % [LA_label.LB1null.h, LA_label.LB1null.rps, LA_label.LB1null.all, LA_label.LB1null.areas] = plotFaceStatistics([LA_label.LB1null.data.Sf]);
% % [LA_label.LB2null.h, LA_label.LB2null.rps, LA_label.LB2null.all, LA_label.LB2null.areas] = plotFaceStatistics([LA_label.LB2null.data.Sf]);
% % [LA_label.mEmerald.h, LA_label.mEmerald.rps, LA_label.mEmerald.all, LA_label.mEmerald.areas] = plotFaceStatistics([LA_label.mEmerald.data.Sf]);
% 
% LA_label.WT.area = getFaceStatistics([LA_label.WT.data.Sf]);
% LA_label.LB1null.area = getFaceStatistics([LA_label.LB1null.data.Sf]);
% LA_label.LB2null.area = getFaceStatistics([LA_label.LB2null.data.Sf]);
% LA_label.mEmerald.area = getFaceStatistics([LA_label.mEmerald.data.Sf]);
% 
% 
% %% Plot QQ-Plot
% figure;
% qqplot(LA_label.WT.area.all.data,LA_label.LB1null.area.all.data)
% xlabel('WT Face Area (\mum^2)')
% ylabel('Lamin B1 null Face Area (\mum^2)')
% title('WT vs LB1null');
% 
% 
% figure;
% qqplot(LA_label.WT.area.all.data,LA_label.LB2null.area.all.data)
% xlabel('WT Face Area (\mum^2)')
% ylabel('Lamin B2 null Face Area (\mum^2)')
% title('WT vs LB2null');
% 
% figure;
% qqplot(LA_label.WT.area.all.data,LA_label.mEmerald.area.all.data)
% xlabel('WT Face Area (\mum^2)')
% ylabel('mEmerald Face Area (\mum^2)')
% title('WT vs mEmerald-LA-fixed');

combined.area.LA.data = [ stats(21).area(1).all.data stats(22).area(1).all.data stats(23).area(2).all.data ];
combined.area.LB1.data = [ stats(21).area(2).all.data stats(25).area(2).all.data stats(27).area(2).all.data ];
combined.area.LB2.data = [ stats(22).area(2).all.data stats(25).area(1).all.data stats(28).area(2).all.data ];
combined.area.LC.data = [ stats(23).area(1).all.data stats(27).area(1).all.data stats(28).area(1).all.data ];
combined.area.all.data = [combined.area.LA.data combined.area.LB1.data combined.area.LB2.data combined.area.LC.data];

lamins.functions.minimalHistogram( {combined.area.LA.data, combined.area.LB1.data, combined.area.LB2.data, combined.area.LC.data } )
xlabel('Face Area (\mum^2)')

[p, tbl, stats] = anova1(combined.area.all.data,[ones(size(combined.area.LA.data)) ones(size(combined.area.LB1.data))*2 ones(size(combined.area.LB2.data))*3 ones(size(combined.area.LC.data))*4])

combined.edgeLength.LA.data = [ stats(21).edgeLength(1).all.data stats(22).edgeLength(1).all.data stats(23).edgeLength(2).all.data ];
combined.edgeLength.LB1.data = [ stats(21).edgeLength(2).all.data stats(25).edgeLength(2).all.data stats(27).edgeLength(2).all.data ];
combined.edgeLength.LB2.data = [ stats(22).edgeLength(2).all.data stats(25).edgeLength(1).all.data stats(28).edgeLength(2).all.data ];
combined.edgeLength.LC.data = [ stats(23).edgeLength(1).all.data stats(27).edgeLength(1).all.data stats(28).edgeLength(1).all.data ];
combined.edgeLength.all.data = [combined.edgeLength.LA.data combined.edgeLength.LB1.data combined.edgeLength.LB2.data combined.edgeLength.LC.data];

lamins.functions.minimalHistogram( {combined.edgeLength.LA.data, combined.edgeLength.LB1.data, combined.edgeLength.LB2.data, combined.edgeLength.LC.data }, (0:40)./sqrt(1000))
xlabel('Edge Length (\mum)')

[p, tbl, stats] = anova1(combined.edgeLength.all.data,[ones(size(combined.edgeLength.LA.data)) ones(size(combined.edgeLength.LB1.data))*2 ones(size(combined.edgeLength.LB2.data))*3 ones(size(combined.edgeLength.LC.data))*4])

combined.edgesPerVertex.LA.data = [ stats(21).edgesPerVertex(1).all.data stats(22).edgesPerVertex(1).all.data stats(23).edgesPerVertex(2).all.data ];
combined.edgesPerVertex.LB1.data = [ stats(21).edgesPerVertex(2).all.data stats(25).edgesPerVertex(2).all.data stats(27).edgesPerVertex(2).all.data ];
combined.edgesPerVertex.LB2.data = [ stats(22).edgesPerVertex(2).all.data stats(25).edgesPerVertex(1).all.data stats(28).edgesPerVertex(2).all.data ];
combined.edgesPerVertex.LC.data = [ stats(23).edgesPerVertex(1).all.data stats(27).edgesPerVertex(1).all.data stats(28).edgesPerVertex(1).all.data ];
combined.edgesPerVertex.all.data = [combined.edgesPerVertex.LA.data combined.edgesPerVertex.LB1.data combined.edgesPerVertex.LB2.data combined.edgesPerVertex.LC.data];

[p, tbl, stats] = anova1(combined.edgesPerVertex.all.data,[ones(size(combined.edgesPerVertex.LA.data)) ones(size(combined.edgesPerVertex.LB1.data))*2 ones(size(combined.edgesPerVertex.LB2.data))*3 ones(size(combined.edgesPerVertex.LC.data))*4])


lamins.functions.minimalHistogram( {combined.edgesPerVertex.LA.data, combined.edgesPerVertex.LB1.data, combined.edgesPerVertex.LB2.data, combined.edgesPerVertex.LC.data } )
xlabel('Edges Per Junction')