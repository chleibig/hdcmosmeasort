function [duplicate_pairs] = CheckForDuplicates(units, sr, ...
                                           t_s, t_jitter, coin_thr,...
                                           sim_thr, interactive)
%CheckForDuplicates tests units pairwisely whether they
%are duplicates or not based on the following criteria:
%
% 1. They have more than coin_thr spikes in common
%
% 2. Their STAs are more similar than sim_thr (similarity is measured
% via their normalized scalar product)
%
% INPUT:
% 
% units: array of structs with at least the fields "time" and "STA"
% sr: sampling rate in kHz
% t_s: maximal overall, pairwise shift of spiketrains to detect coincident
%      spikes in ms
% t_jitter: two spikes are considered as coincident if they are maximally
%           t_s + t_jitter ms apart from each other 
% coin_thr: minimum fraction of coincident spikes
% sim_thr: minimum similarity of average waveforms (STA)
% 
% OUTPUT:
%
% duplicate_pairs: N_duplicates x 2 matrix containing the indices that 
%                  were found to be duplicates; the first index is the
%                  one that contributes with the higher fraction of coinci-
%                  dent spikes

% author: Christian Leibig, 19.03.13

%procedure: we check all pairs for the most evident duplicate pair, remove
%           that one and then do the same for the rest and so on until
%           we do not find any more duplicate
%           NOTE: the way it is currently implemented is ordered in the
%           following fashion: first look for the most coincident pair and
%           then look for the most similar STA, is the other way
%           around more preferable? Consider signal-to-noise ratio optimal
%           scalar-product to assess similarity

if interactive
   figure; 
   %this is for debug purposes and to assess suitable threshold parameters
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
   subplot(1,2,1);imagesc(coin_frac);colorbar;
   title('coincident spikes');
   xlabel('unit index');
   ylabel('unit index');
   axis square;
   sim = diag(ones(N,1));
   for i=1:N
       for j=i+1:N
           STAi = units(i).STA(:);
           STAj = units(j).STA(:);
           sim(i,j) = abs(STAi'*STAj / sqrt((STAi'*STAi) * (STAj'*STAj)));
           sim(j,i) = sim(i,j);
       end
   end
   subplot(1,2,2);imagesc(sim);colorbar;
   title('similarities of STAs');
   xlabel('unit index');
   ylabel('unit index');
   axis square;
   coin_thr = input('Please specify the minimum fraction of coincident spikes for duplicates:');
   sim_thr = input('Please specify the minimum similarity of STAs for duplicates:');
   
   clear N sim coin_frac
end
    

units_tmp = units; %units_tmp will be changed
%store the original indices away:
for i = 1:length(units_tmp); units_tmp(i).ID = i; end

duplicate_pairs = [];

continue_to_check = true;

while continue_to_check

    continue_to_check = false;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Coincident spikes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    N = length(units_tmp);
    coin_frac = zeros(N,N);
    for i=1:N
        for j=i+1:N
            ti = units_tmp(i).time;
            tj = units_tmp(j).time;
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
        continue_to_check = false;
    end
    
    if ~isempty(I)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STA similarities
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        STAi = units_tmp(I).STA(:);
        STAj = units_tmp(J).STA(:);
        sim = abs(STAi'*STAj / sqrt((STAi'*STAi) * (STAj'*STAj)));
        if sim >= sim_thr
            %A duplicate pair is found
            duplicate_pairs = [duplicate_pairs;...
                              units_tmp(I).ID units_tmp(J).ID];
            units_tmp(I) = [];%exclude that from further checks
            continue_to_check = true;
        else
            %this was not a duplicate, but we only stop if both criteria
            %are false
            units_tmp(I) = [];%exclude that from further checks
            continue_to_check = true;
        end
    
    end    
end

end