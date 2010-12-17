function [groupNeigh]=groupNeighResults(groupNeigh,neigh)
[grpBin,grpFrame]=size(groupNeigh);
[sglBin,sglFrame]=size(     neigh);

maxBin  =max(grpBin  ,sglBin);
maxFrame=max(grpFrame,sglFrame);

% c2cDis c2cDD dirDis trvDis c2edMean c2edSTD

for bin=1:maxBin
    for frame=1:maxFrame
        if bin>sglBin || frame>sglFrame || isempty(neigh(bin,frame))
            neigh=initNeigh(neigh,bin,frame);
        end
        if bin>grpBin || frame>grpFrame || isempty(groupNeigh(bin,frame))
            groupNeigh=initNeigh(groupNeigh,bin,frame);
        end        
        
        grpNumDt=size(groupNeigh(bin,frame).c2cDis,2);
        newNumDt=size(     neigh(bin,frame).c2cDis,2);
        
        groupNeigh(bin,frame).c2cDis  = vertcat(padarray(groupNeigh(bin,frame).c2cDis , [0,max(grpNumDt,newNumDt)-grpNumDt], NaN  ,'post'), padarray(     neigh(bin,frame).c2cDis , [0,max(grpNumDt,newNumDt)-newNumDt], NaN  ,'post'));
        groupNeigh(bin,frame).c2cDD   = vertcat(padarray(groupNeigh(bin,frame).c2cDD  , [0,max(grpNumDt,newNumDt)-grpNumDt], NaN  ,'post'), padarray(     neigh(bin,frame).c2cDD  , [0,max(grpNumDt,newNumDt)-newNumDt], NaN  ,'post'));
        groupNeigh(bin,frame).dirDis  = vertcat(padarray(groupNeigh(bin,frame).dirDis , [0,max(grpNumDt,newNumDt)-grpNumDt], NaN  ,'post'), padarray(     neigh(bin,frame).dirDis , [0,max(grpNumDt,newNumDt)-newNumDt], NaN  ,'post'));
        groupNeigh(bin,frame).trvDis  = vertcat(padarray(groupNeigh(bin,frame).trvDis , [0,max(grpNumDt,newNumDt)-grpNumDt], NaN  ,'post'), padarray(     neigh(bin,frame).trvDis , [0,max(grpNumDt,newNumDt)-newNumDt], NaN  ,'post'));
        % the in-sheet positions:
        groupNeigh(bin,frame).c2edMean= vertcat(groupNeigh(bin,frame).c2edMean,neigh(bin,frame).c2edMean);
        groupNeigh(bin,frame).c2edSTD = vertcat(groupNeigh(bin,frame).c2edSTD ,neigh(bin,frame).c2edSTD);
        
        
        
%         % also store the full vectors:
%         grpR=groupNeigh(bin,frame).R;
%         sglR=neigh(bin,frame).R;
%         
%         grpc2cBin=length(grpR);
%         sglc2cBin=length(sglR);        
%         
%         maxc2cBin=max(grpc2cBin,sglc2cBin);
%         for c2cBin=1:maxc2cBin
%             if c2cBin>sglc2cBin || isempty(sglR(c2cBin))
%                 sglR=initR(sglR,c2cBin);
%             end
%             if c2cBin>grpc2cBin || isempty(grpR(c2cBin))
%                 grpR=initR(grpR,c2cBin);
%             end           
%             
%             grpR(c2cBin).pts2ptR = vertcat(grpR(c2cBin).pts2ptR,sglR(c2cBin).pts2ptR);
%             grpR(c2cBin).xVelPts = vertcat(grpR(c2cBin).xVelPts,sglR(c2cBin).xVelPts);
%             grpR(c2cBin).yVelPts = vertcat(grpR(c2cBin).yVelPts,sglR(c2cBin).yVelPts);
%             grpR(c2cBin).xVelCtr = vertcat(grpR(c2cBin).xVelCtr,sglR(c2cBin).xVelCtr);
%             grpR(c2cBin).yVelCtr = vertcat(grpR(c2cBin).yVelCtr,sglR(c2cBin).yVelCtr);
%             grpR(c2cBin).cosda   = vertcat(grpR(c2cBin).cosda  ,sglR(c2cBin).cosda);
%             % increase the total number of correlation pairs:
%             grpR(c2cBin).Ntot    = vertcat(grpR(c2cBin).Ntot,sglR(c2cBin).Ntot);
%         end
%         groupNeigh(bin,frame).R=grpR;
    end
end

function newNeigh=initNeigh(newNeigh,newBin,newFrame)
newNeigh(newBin,newFrame).c2cDis          = [];
newNeigh(newBin,newFrame).c2cDD           = [];
newNeigh(newBin,newFrame).dirDis        = [];
newNeigh(newBin,newFrame).trvDis       = [];
% the in-sheet positions:
newNeigh(newBin,newFrame).c2edMean       = [];
newNeigh(newBin,newFrame).c2edSTD        = [];



% function newR=initR(newR,newRBin)
% newR(newRBin).pts2ptR = [];
% newR(newRBin).xVelPts = [];
% newR(newRBin).yVelPts = [];
% newR(newRBin).xVelCtr = [];
% newR(newRBin).yVelCtr = [];
% newR(newRBin).cosda   = [];
% % increase the total number of correlation pairs:
% newR(newRBin).Ntot    = [];