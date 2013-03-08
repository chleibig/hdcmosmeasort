function [I,J] = get_duplicate(X,sr,t_s,t_jitter, coin_thr)
%GET_DUPLICATE returns the duplicate index pair for the most evident
% duplicate; currently "duplicateness" is assessed via the fraction
% of coincident spikes. I is the component that contains less spikes than
% J.
%
% Input arguments:
%
% X: (component x samples) array
% sr: sampling rate in kHz
% t_s: maximal global shift of spike trains against each other in ms
% t_jitter: maximal additional distance between two spikes to be counted as
%           coincident in ms
% NOTE: coincident spikes are maximally t_s + t_jitter ms apart
% coin_thr: minimal fraction of coincident spikes for a pair to be
%           considered as duplicate

% Author: Christian Leibig, 08.03.13

[units] = SpikeTimeIdentification(X, sr, 0);
N = length(units);
coin_frac = zeros(N,N);
for i=1:N
    for j=i+1:N
        ti = units(i).time;
        tj = units(j).time;
        [n_coin] = count_coincident_spks(ti,tj,sr,t_s,t_jitter);
        coin_frac(i,j) = n_coin/length(ti);
        coin_frac(j,i) = n_coin/length(tj);
    end
end

[val_max, arg_max] = max(coin_frac(:));
if val_max >= coin_thr
    %I corresponds to the component to which coin_frac
    %is normalized for the maximum value
    [I,J] = ind2sub(size(coin_frac),arg_max);
else
    I = [];J = [];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Incorporation of mixing filter / template similarity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% due to linear filtering ambiguity in case of convolutive mixing it is not
% straight forward to calculate the similarity based on the mixing
% matrices, because the maximum signal energy might be at different lags tau.
% A cleaner way would be, to go back to the raw data and get the
% templates with the spike times; however also there proper alignment has
% to be guaranteed

%test with incorporating only the first (mainly dominated by fastICA)
%mixing matrix gave only low similarities (~ 0.1 for coin_frac = 1!)
%     Mij = abs(A_tau(:,I,1)'*A_tau(:,J,1) / ...
%         (sqrt(A_tau(:,I,1)'*A_tau(:,I,1)) *
%         sqrt(A_tau(:,J,1)'*A_tau(:,J,1))));

end



