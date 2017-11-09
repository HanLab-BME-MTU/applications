%% Setup
import lamins.functions.*;
import lamins.classes.*;

path = '/project/biophysics/jaqaman_lab/lamins/2015/20150602/MEFLB1-LACLB12/MEFLB1-LACLB12-006_Reconstructed';
MD = MovieData.load([path filesep 'MEFLB1-LACLB12-006_Reconstructed.mat']);
% cd(path);
skeletons = load([MD.outputDirectory_ filesep 'skeletons_2015_06_10.mat']);

L = LaminsData(MD);
images = L.getImages;
ctz = lamins.util.TZtoCTZ_LinearInd(MD,skeletons.tz);
I = images(ctz(1));

%% Original Image

name = '00_image';
h = figure;
imshow(I);
set(h,'PaperPosition',[0 0 1 1]);
print(h,name,'-dpng','-r400');
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-dpng','-r400');
close(h);

%% Green Image

name = '00_image_green';
h = figure;
rgb = zeros(1024,1024,3);
rgb(:,:,2) = double(I);
imshow(rgb);
set(h,'PaperPosition',[0 0 1 1]);
print(h,name,'-dpng','-r400');
set(h,'PaperPosition',[0 0 6 6]);
print(h,[name '_six'],'-dpng','-r400');
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-dpng','-r400');
close(h);

%% Flattened Image

name = '01_flattened';
h = figure;
F = I.flattenIntensity;
imshow(;
set(h,'PaperPosition',[0 0 1 1]);
print(h,name,'-dpng','-r400');
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-dpng','-r400');
close(h);


%% Steerable response

name = '02_steerable';
h = figure;
imshow(I.steerable.res,[]);
set(h,'PaperPosition',[0 0 1 1]);
print(h,name,'-dpng','-r400');
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-dpng','-r400');
close(h);


%% Steerable nms

name = '03_steerable';
h = figure;
imshow(I.steerable.nms,[]);
set(h,'PaperPosition',[0 0 1 1]);
print(h,name,'-dpng','-r400');
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-dpng','-r400');
close(h);

%% Steerable nmsSkel

name = '04_nmsSkel';
h = figure;
imshow(I.nmsSkel,[]);
set(h,'PaperPosition',[0 0 1 1]);
print(h,name,'-dpng','-r400');
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-dpng','-r400');
close(h);


%% Thresh Skel

name = '05_threshSkel';
h = figure;
imshow(I.threshSkel,[]);
set(h,'PaperPosition',[0 0 1 1]);
print(h,name,'-dpng','-r400');
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-dpng','-r400');
close(h);


%% Mask

name = '06_mask';
h = figure;
imshow(I.mask,[]);
set(h,'PaperPosition',[0 0 1 1]);
print(h,name,'-dpng','-r400');
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-dpng','-r400');
close(h);


%% Mask on Image

name = '07_mask_image';
h = figure;
X = I.mask.* im2double(I.adjusted);
X = X + ~I.mask*0.5;
imshow(X);
set(h,'PaperPosition',[0 0 1 1]);
print(h,name,'-dpng','-r400');
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-dpng','-r400');
close(h);


%% Extended Skel

name = '08_extendedSkel';
h = figure;
imshow(I.extendedSkel,[]);
set(h,'PaperPosition',[0 0 1 1]);
print(h,name,'-dpng','-r400');
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-dpng','-r400');
close(h);

%% Skeleton

S = I.skeleton;

name = '09_skeleton';
h = figure;
imshow(S);
set(h,'PaperPosition',[0 0 1 1]);
print(h,name,'-dpng','-r400');
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-dpng','-r400');
close(h);

%% Skeleton convertShortEdgesToVertices

name = '10_convertShortEdgesToVertices';
S.convertShortEdgesToVertices(2);
h = figure;
imshow(S);
set(h,'PaperPosition',[0 0 1 1]);
print(h,name,'-dpng','-r400');
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-dpng','-r400');
close(h);

%% Skeleton reduceVerticestoPoints

name = '11_reduceVerticesToPoints';
S.reduceVerticesToPoints;
h = figure;
imshow(S);
set(h,'PaperPosition',[0 0 1 1]);
print(h,name,'-dpng','-r400');
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-dpng','-r400');
close(h);

%% S drawEdgesAslines

name = '12_drawEdgesAsLines';
h = figure;
imshow(I);
% Sb.drawEdgesAsLines([],'g');
S.drawEdgesAsLines([],'m');
set(h,'PaperPosition',[0 0 6 6]);
print(h,name,'-depsc2','-r400');
set(h,'PaperPosition',[0 0 1 1]);
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-depsc2','-r400');
close(h);
%% Audit by Mask
% Auditing 1
name = '13_auditByMask';
Sb = S.copy;

                % New on June 3rd, 2015
                % Audit using mask, proportion on flattened intensity, and
                % do another round of score optimization including zero
                S.auditEdgesByMask(I);
                
h = figure;
imshow(I);
Sb.drawEdgesAsLines([],'g');
S.drawEdgesAsLines([],'m');
set(h,'PaperPosition',[0 0 6 6]);
print(h,name,'-depsc2','-r400');
set(h,'PaperPosition',[0 0 1 1]);
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-depsc2','-r400');
close(h);

%% Skeleton removeEdgesWithoutFace
% Auditing 2
name = '14_removeEdgesWithoutFace';
Sb = S.copy;
            % Cleans up the edges and faces of the skeleton
            % 1. Removes edges that have no faces
            % 2. Removes faces that have no edges
            % Future considerations: Remove small faces
            [~,EF] = S.faceEdges;
            e = find(cellfun(@length,EF) == 0);
            e = [ e ; (length(EF)+1:S.edges.NumObjects)'];
            S.deleteEdges(e);
            FE = S.faceEdges;
            f = find(cellfun(@length,FE) == 0);
            S.deleteFaces(f);
h = figure;
imshow(I);
Sb.drawEdgesAsLines([],'g');
S.drawEdgesAsLines([],'m');
set(h,'PaperPosition',[0 0 6 6]);
print(h,name,'-depsc2','-r400');
set(h,'PaperPosition',[0 0 1 1]);
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-depsc2','-r400');
close(h);



%% Score and delete 1
% Auditing 3
name = '15_scoreAndDeleteFirst';
Sb = S.copy;
                score = lamins.functions.scoreEdges(S,I.flattenIntensity);
                S.deleteEdges(score < 0);
                
h = figure;
imshow(I);
Sb.drawEdgesAsLines([],'g');
S.drawEdgesAsLines([],'m');
set(h,'PaperPosition',[0 0 6 6]);
print(h,name,'-depsc2','-r400');
set(h,'PaperPosition',[0 0 1 1]);
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-depsc2','-r400');
close(h);

%% Audit according to intensity variation
% Auditing 4 and 5
name = '16_auditAccordingToIntensityVariation';
Sb = S.copy;
                % Use intensity variation along the edge
                S2 = S.copy;
                thresh = I.maskThresh(double(I));
                S2.auditEdges(double(I),[],thresh,thresh/2);
                
h = figure;
imshow(I);
Sb.drawEdgesAsLines([],'g');
S2.drawEdgesAsLines([],'m');
set(h,'PaperPosition',[0 0 6 6]);
drawnow;
print(h,name,'-depsc2','-r400');
set(h,'PaperPosition',[0 0 1 1]);
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-depsc2','-r400');
close(h);



%% Audit edges by thresolded intensity
% Auditing 6
Sb = S2.copy;
                S3 = S2.copy;
                S3.auditEdgesByThresholdedIntensity(I);

name = '17_auditByThresholdedIntensity';


h = figure;
imshow(I);
Sb.drawEdgesAsLines([],'g');
S3.drawEdgesAsLines([],'m');
set(h,'PaperPosition',[0 0 6 6]);
print(h,name,'-depsc2','-r400');
drawnow;
set(h,'PaperPosition',[0 0 1 1]);
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-depsc2','-r400');
close(h);

%% Score and delete 2
% Auditing 7

name = '18_scoreAndDeleteSecond';
Sb = S3.copy;

                
                score2 = lamins.functions.scoreEdges(S3,I.flattenIntensity);
                S3.deleteEdges(score2 <= 0);
                
h = figure;
imshow(I);
Sb.drawEdgesAsLines([],'g');
S3.drawEdgesAsLines([],'m');
set(h,'PaperPosition',[0 0 6 6]);
drawnow;
print(h,name,'-depsc2','-r400');
set(h,'PaperPosition',[0 0 1 1]);
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-depsc2','-r400');
close(h);

%% Final

name = '19_final';
Sb = S3.copy;

                
                
h = figure;
imshow(I);
S3.drawEdgesAsLines([],'m');
set(h,'PaperPosition',[0 0 6 6]);
print(h,name,'-depsc2','-r400');
drawnow;
set(h,'PaperPosition',[0 0 1 1]);
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-depsc2','-r400');
close(h);

%% Scale

name = 'scale';

Z = zeros(1024);
Z(453:652,306:505) = 1;
Z(564,340:(340+159)) = 0;

h = figure;
imshow(Z);
set(h,'PaperPosition',[0 0 1 1]);
print(h,name,'-dpng','-r400');
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-dpng','-r400');
close(h);

%% Scale EPS

name = 'scale';

Z = zeros(1024);
Z(453:652,306:505) = 1;
Z(564,340:(340+159)) = 0;

h = figure;
imshow(Z);
set(h,'PaperPosition',[0 0 6 6]);
print(h,name,'-depsc2','-r400');
drawnow;
set(h,'PaperPosition',[0 0 1 1]);
xlim([306 505]);
ylim([453.5 652.5]);
print(h,[name '_zoom'],'-depsc2','-r400');
close(h);

