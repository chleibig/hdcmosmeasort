function [n_coin] = count_coincident_spks(t1,t2,sr,t_s,t_jitter)
%count_coincident_spks counts coincident spikes of spike trains t1 and t2
%
% Input:
%
% t1, t2: vectors of timestamps in ms
% sr: sampling rate in kHz
% t_s: maximal (global) shift of spike trains against each other [ms]
% t_jitter: maximal distance between pairs of spikes in ms for the shifted
%           spiketrain
% maximally allowed distance between spikes to be counted as coincident 
% events is hence given by dt = t_s + t_jitter

% author: Christian Leibig, 07.03.2013

% Ref.: Felix Franke, 2012, PhD Thesis, chapter 7.1


%shift spike train by maximally t_s:
edges = min(t1(1),t2(1)):1/sr:max(t1(end),t2(end));
t1_binned = histc(t1,edges);
t2_binned = histc(t2,edges);
maxlags = ceil(t_s*sr);%in frames
c = xcorr(t1_binned,t2_binned,maxlags,'coeff');
[max_value,max_shift] = max(c);
t2_shift = (max_shift - (maxlags + 1) )/sr; % in ms
t2 = t2 + t2_shift;

%count coincident spikes:
n_coin = 0;
for k=1:length(t1)
    if any((t2 - t_jitter <= t1(k)) & ...
            (t1(k) <= t2 + t_jitter))
        n_coin = n_coin + 1;
    end
end
    
    
    
end