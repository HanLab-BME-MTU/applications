% clear;
% 
% M = imread('C:\Documents and Settings\godina\Desktop\Skeletons\Movies\Series4\nakedSkeletons\Four0100.tif');
% mat = M; 
% 
% [ConnectivityMatrix,ConnectivityMatrixRussianForm,PositionVectices] = GetConnectivityPointbyPoint(mat);
% 
% x = PositionVectices(:,1);
% y = PositionVectices(:,2);
% tic
%     [dSP,sp]=grShortPath(ConnectivityMatrixRussianForm,1,200);
% %     round(size(ConnectivityMatrix,1)));
% toc
% 
% A = zeros(size(mat,1),size(mat,2),3);
% A(:,:,2) = mat;
% for it = 1:length(sp)
% A(x(sp(it)),y(sp(it)),1) = 1;
% end
% seeImage(A);
% %%
% % 
% % i1 = 8;
% % i2 = 4;
% % 
% % lo = 1;
% % Ctemp = C;
% % while Ctemp(i1,i2)==0
% %     lo = lo + 1
% %     tic
% %     Ctemp = C.^lo;
% %     seeImage(Ctemp)
% %     toc
% %     [lo,length(find(Ctemp))]
% % end
% % %%
% % C = [0 0 0 0 0 0;
% %      1 0 0 0 0 0;
% %      0 1 0 0 0 0;
% %      0 0 0 0 0 0;
% %      1 0 0 0 0 0;
% %      0 0 0 0 0 0];
% % Con = C + C'
% % 
% % X = 2;
% % Y = 5;
% % 
h = waitbar(.3,'aal')
h = waitbar(.5,h,'aal  2')