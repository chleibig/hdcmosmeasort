function [n_coin, syncSpks1, syncSpks2] = spiketrainalignment(t1,t2,sr,t_s,t_jitter)
%spiketrainalignment counts coincident spikes of spike trains t1 and t2
%and tracks which spikes in t1 and t2 participate in synchronous events
% 
% Input
% =====
%
% t1, t2: vectors of timestamps in ms
% sr: sampling rate in kHz
% t_s: maximal (global) shift of spike trains against each other [ms]
% t_jitter: maximal distance between pairs of spikes in ms for the shifted
%           spiketrain
% maximally allowed distance between spikes to be counted as coincident 
% events is hence given by dt = t_s + t_jitter
%
%
% Output
% ======
%
% n_coin - the number of coincident spikes
% syncSpks1, *2 - indices into t1 and t2 indicating respective spikes which
%                 participate in coincident events

% author: christian.leibig@g-node.org, 07.03.2013, updated 28.08.13
% 

% Ref.: Felix Franke, 2012, PhD Thesis, chapter 7.1


if isempty(t1) || isempty(t2); n_coin = 0; return ;end

%shift spike train by maximally t_s:
edges = min(t1(1),t2(1)):1/sr:max(t1(end),t2(end));
t1_binned = histc(t1,edges);
t2_binned = histc(t2,edges);
maxlags = ceil(t_s*sr);%in frames
c = xcorr(t1_binned,t2_binned,maxlags,'coeff');
[unused,max_shift] = max(c);
t2_shift = (max_shift - (maxlags + 1) )/sr; % in ms
t2 = t2 + t2_shift;


%count coincident spikes and track which ones participate in coincident
%events:
n_coin = 0;
syncSpks1 = [];
syncSpks2 = [];
for k=1:length(t1)
    t2Lower = find((t2 - t_jitter) <= t1(k), 1,'last');
    t2Upper = find(t1(k) <= (t2 + t_jitter), 1,'first');
    if ~isempty(t2Lower) && ~isempty(t2Upper) && (t2Upper - t2Lower == 0)
        n_coin = n_coin + 1;
        syncSpks1 = [syncSpks1;k];
        syncSpks2 = [syncSpks2;t2Lower];      
    end
end

    
end