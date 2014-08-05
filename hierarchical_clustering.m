function [T] = hierarchical_clustering(SM,min_sim,varargin)
%Wrapper function for hierarchical clustering of Statistics Toolbox

%Input:

%SM upper triangular similarity matrix

% author: Christian Leibig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotting = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optional arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (rem(length(varargin),2)==1)
  error('Optional parameters should always go by pairs');
else
  for i=1:2:(length(varargin)-1)
    if ~ischar (varargin{i}),
      error (['Unknown type of optional parameter name (parameter' ...
	      ' names must be strings).']);
    end
    % change the value of parameter
    switch varargin{i}
     case 'plotting'
      plotting = (varargin{i+1});     
     otherwise
      % Hmmm, something wrong with the parameter string
      error(['Unrecognized parameter: ''' varargin{i} '''']);
     end;
  end;
end

% Check if there is anything to cluster at all:
if length(SM) < 2
    fprintf('Warning: Clustering skipped because less than two input identities are provided.\n');
    T = 1;
    return
end

%Conversion of similarity matrix SM to output format of pdist
% Y = reshape(SM', [1 size(SM,1)*size(SM,2)]);
% Y = nonzeros(Y);
N = length(SM);
Y = zeros([1 N*(N-1)/2]);
for row = 1:N-1
    for col = row + 1
        %offset analytically from difference of geometric series:
        offset = N*(row - 1) + row/2 - row^2/2;
        %more readable: offset = N*(N-1)/2 - sum(1:(N-row))
        Y(offset+1:offset+N-row) = SM(row,col:end);
    end
end

%Conversion from similarities to distances:
Y = 1 - Y;
max_dist = 1 - min_sim;
% Z = linkage(Y','single');
Z = linkage(Y,'single');
if plotting
    figure;
    subplot(2,2,1);imagesc(SM);axis square; title('crosstalk');cb = colorbar;
    set(cb,'CLim',[0,1]);
    [I,J] = find(abs(SM) > min_sim);hold on;scatter(J,I);
    subplot(2,2,2);hist(SM(find(triu(SM,1))),sqrt(nnz(triu(SM,1))));
    title('hist(crosstalk)');
    xlabel('crosstalk');ylabel('counts');
    subplot(2,2,3:4);dendrogram(Z);title('dendrogram based on distances');
    ylim([0,1]);
    hold on; plot([0 size(SM,1)+1],[max_dist max_dist],'r');
%     min_sim_tmp = input(strcat('Current min. crosstalk is ',num2str(min_sim),...
%     '.\nIf desired, please enter different threshold: '));
%     if ~isempty(min_sim_tmp);max_dist = 1 - min_sim_tmp;end
%     clear min_sim_tmp
end

T = cluster(Z,'Cutoff',max_dist,'Criterion','distance');

end
