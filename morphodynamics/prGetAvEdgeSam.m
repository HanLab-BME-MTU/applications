function [x_normal_out_av, y_normal_out_av, x_av_pos_normal, y_av_pos_normal] = prGetAvEdgeSam(matrix, varargin)
% PRGETAVEDGE averages vectors assiciated along an edge
% 
%             Given a spline and values along that spline, this function
%             devides the spline in a number of segments and calculates the
%             average of this vectors in each segment.
%             It is assumed that the spline parameter starts with 1 and
%             ends with max_spline_par. If there is no mean value found for
%             a segment the values are set to NaN.
%             
%
% SYNOPSIS    [x_av, y_av, x_n_av, y_n_av, pos] = prGetAvEdge(sp_x, sp_y, max_spline_par, p_vec, x_vec, y_vec, varargin)
%
% INPUT       matrix                     : protrusion or normal matrix (x_pos, y_pos, x_comp, y_comp)
% 
% OUTPUT      x_normal_out_av            : x-comp of averaged vector
%             y_normal_out_av            : y-comp of averaged vector
%             x_av_pos_normal            : x-comp of averaged vector position
%             y_av_pos_normal            : x-comp of averaged vector position
%                           
% DEPENDENCES   prGetAvEdgeSam uses { 
%                                 } 
%               prGetAvEdgeSam is used by { imEdgeTracker
%                                 } 
%
% Shann-Ching Chen 03/20/2008

%         [x_normal_out_av, y_normal_out_av, x_av_pos_normal, y_av_pos_normal, m_pos]=...
%             prGetAvEdge(edge_sp_array_x(time), edge_sp_array_y(time), edge_sp_array_x(time).knots(end),...
%             p_n, normal_matrix{time}(:,3), normal_matrix{time}(:,4),
%             'nr_sect', SEG_NR);


%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=length(varargin);

for i=1:l
   if strcmp(varargin(i),'nr_sect')
      NR_SECT=varargin{i+1};    
   end
end

if ~exist('NR_SECT','var')
   NR_SECT=10;
end
%%%%%%%%%%%%%%%%%%% End parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the spline parameter interval based
% on given number of sections
nfold = NR_SECT;
pixelNum = size(matrix,1);
binsize = floor(pixelNum/nfold);
%halfbin = round(binsize/2);
%start_idx = 1 - halfbin + pixelNum;
start_idx = 0;
idx{1} = mod(start_idx:(start_idx+binsize-1), pixelNum)+1;
start_idx = mod(start_idx+binsize, pixelNum)+1;
for w=2:(nfold-mod(pixelNum, nfold))
    idx{w} = start_idx:(start_idx+binsize-1);
    start_idx  = start_idx+binsize;
end
for w=(nfold+1-mod(pixelNum, nfold)):nfold
    idx{w} = start_idx:(start_idx+binsize);
    start_idx  = start_idx+binsize+1;
end

% average vectors along the edge
for i=1:NR_SECT
    x_av_pos_normal(i) = 0;
    y_av_pos_normal(i) = 0;    
    x_normal_out_av(i) = 0;
    y_normal_out_av(i) = 0;

    for j=1:length(idx{i})
        x_av_pos_normal(i) = x_av_pos_normal(i)+matrix(idx{i}(j),1);
        y_av_pos_normal(i) = y_av_pos_normal(i)+matrix(idx{i}(j),2);        
        x_normal_out_av(i) = x_normal_out_av(i)+matrix(idx{i}(j),3);
        y_normal_out_av(i) = y_normal_out_av(i)+matrix(idx{i}(j),4);
    end
    
    if length(idx{i}) > 0 
        % averaged values
        x_av_pos_normal(i)=x_av_pos_normal(i)/length(idx{i});
        y_av_pos_normal(i)=y_av_pos_normal(i)/length(idx{i});      
        x_normal_out_av(i)=x_normal_out_av(i)/length(idx{i});
        y_normal_out_av(i)=y_normal_out_av(i)/length(idx{i});
    else
        x_av_pos_normal(i)=NaN;
        y_av_pos_normal(i)=NaN;       
        x_normal_out_av(i)=NaN;
        y_normal_out_av(i)=NaN;   
    end
end
