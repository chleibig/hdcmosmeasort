function [numDistinct,sizes,members] = checkforredundancy(units, sr, ...
                                           t_s, t_jitter, coin_thr,...
                                           sim_thr, show, interactive,varargin)
% [numDistinct,sizes,members] = checkforredundancy(units, sr, ...
%                                    t_s, t_jitter, coin_thr,...
%                                      sim_thr, show, interactive,varargin)
% groups units by redundancy of any order (pairwise or higher) based on the 
% following criteria
%
% 1. They have (pairwise) more than coin_thr spikes in common
%
% 2. Their STAs are more similar than sim_thr (similarity is measured
% via their pairwise normalized scalar product)
%
% Both pairwise criteria matrices from 1. and 2. are thresholded by
% coin_thr and sim_thr respectively, then logically anded into a single
% matrix. The resulting binary matrix is then decomposed into numDistinct
% connected components. Each connected component corresponds
% to a distinct underlying unit with which one or more members are
% associated
% 
%
% INPUT
% =====
%
% units: array of structs with at least the fields "time" and "STA"
% sr: sampling rate in kHz
% t_s: maximal overall, pairwise shift of spiketrains to detect coincident
%      spikes in ms
% t_jitter: two spikes are considered as coincident if they are maximally
%           t_s + t_jitter ms apart from each other 
% coin_thr: minimum fraction of coincident spikes, both participants in
%           pairwise comparison have to exceed coin_thr
% sim_thr: minimum similarity of average waveforms (STA)
%
% Optional Input 
% ==============
%
% 'unitIDs',unitIDs - vector of unit indices
% 
% OUTPUT
% ======
%
% numDistinct - number of distinct underlying components
% sizes - vector, sizes(i) specifies the number of members of connected
%         distinct unit i
% members - cell array, members{i} contains the indices of the original
%           units associated with the ith distinct unit

% christian.leibig@g-node.org, 10.10.2013


N = length(units);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

unitIDs = 1:N;

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
        case 'unitIDs'
            unitIDs = (varargin{i+1});
        otherwise
            % Hmmm, something wrong with the parameter string
            error(['Unrecognized parameter: ''' varargin{i} '''']);
    end;
  end;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computation of N(N - 1)/2 pairwise criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
% Coincident spikes %
%%%%%%%%%%%%%%%%%%%%%

coin_frac = zeros(N,N);
for i=1:N
    for j=i+1:N
        ti = units(i).time;
        tj = units(j).time;
        [n_coin,unused1,unused2] = spiketrainalignment(ti,tj,sr,t_s,t_jitter);
        coin_frac(i,j) = n_coin/length(ti);
        coin_frac(j,i) = n_coin/length(tj);
    end
end

%%%%%%%%%%%%%%%%%%%%
% STA similarities %
%%%%%%%%%%%%%%%%%%%%

sim = zeros(N,N);
for i=1:N
    for j=i+1:N
        %temporal alignment on peaks:
        STAi = reshape(units(i).STA,...
         [size(units(i).STA,1)*size(units(i).STA,2) size(units(i).STA,3)]);
        STAj = reshape(units(j).STA,...
         [size(units(j).STA,1)*size(units(j).STA,2) size(units(j).STA,3)]);
        [i_v_p,i_p] = max(max(abs(STAi),[],1));
        [j_v_p, j_p] = max(max(abs(STAj),[],1));
        shift = abs(i_p - j_p);
        if i_p >= j_p
            STAi = STAi(:,1+shift:end);
            STAj = STAj(:,1:end-shift);
        else
            STAj = STAj(:,1+shift:end);
            STAi = STAi(:,1:end-shift);
        end
        STAi = STAi(:);
        STAj = STAj(:);
        %calculate similarity:
        sim(i,j) = STAi'*STAj / sqrt((STAi'*STAi) * (STAj'*STAj));
        sim(j,i) = sim(i,j);
    end
end


if show
    fig1 = figure;
    subplot(3,2,1);imagesc(coin_frac);colorbar;
    title('Coincident spikes');
    xlabel('unit index');
    ylabel('unit index');
    axis square;
    set(gca,'XTick',1:N);
    set(gca,'XTickLabel',unitIDs);
    set(gca,'YTick',1:N);
    set(gca,'YTickLabel',unitIDs);

    subplot(3,2,2);hist(coin_frac(:));
    xlabel('fraction of coincident spikes');
    ylabel('counts');

    subplot(3,2,3);imagesc(sim);colorbar;
    title('Similarities of STAs');
    xlabel('unit index');
    ylabel('unit index');
    axis square;
    set(gca,'XTick',1:N);
    set(gca,'XTickLabel',unitIDs);
    set(gca,'YTick',1:N);
    set(gca,'YTickLabel',unitIDs);

    subplot(3,2,4);hist(sim(:));
    xlabel('c_{ij}');
    ylabel('counts');
end

if interactive
   coin_thr = input('Please specify the minimum fraction of coincident spikes for duplicates:');
   sim_thr = input('Please specify the minimum similarity of STAs for duplicates:');
end
    
%coin_frac is generally asymmetric - ensure that for every entry as well
%the entry at the transposed position is above threshold.
combinedCriteria = ((coin_frac >= coin_thr) & (coin_frac >= coin_thr)')...
                   & (sim >= sim_thr);

if show
    figure(fig1);
    subplot(3,2,5);imagesc(combinedCriteria);colorbar;
    title('Combined criteria');
    xlabel('unit index');
    ylabel('unit index');
    axis square;
    set(gca,'XTick',1:N);
    set(gca,'XTickLabel',unitIDs);
    set(gca,'YTick',1:N);
    set(gca,'YTickLabel',unitIDs);
end

[numDistinct,sizes,members,unused] = networkComponents(combinedCriteria);

%Map member indices to unitIDs.
members = cellfun(@(x) unitIDs(x),members,'UniformOutput',false);

end