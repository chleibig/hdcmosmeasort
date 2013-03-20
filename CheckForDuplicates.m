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
%


   
N = length(units);

%first compute all the N(N - 1)/2 pairwise criteria:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coincident spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STA similarities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sim = diag(ones(N,1));
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
        sim(i,j) = abs(STAi'*STAj / sqrt((STAi'*STAi) * (STAj'*STAj)));
        sim(j,i) = sim(i,j);
    end
end


if interactive
   figure; 
   %this is for debug purposes and to assess suitable threshold parameters
   subplot(1,2,1);imagesc(coin_frac);colorbar;
   title('coincident spikes');
   xlabel('unit index');
   ylabel('unit index');
   axis square;
   subplot(1,2,2);imagesc(sim);colorbar;
   title('similarities of STAs');
   xlabel('unit index');
   ylabel('unit index');
   axis square;
   coin_thr = input('Please specify the minimum fraction of coincident spikes for duplicates:');
   sim_thr = input('Please specify the minimum similarity of STAs for duplicates:');
end
    
%*_tmp criteria matrices will be changed in place to keep track of which
%units were already considered as duplicates and should be excluded from
%further comparisons:
coin_frac_tmp = coin_frac; 
sim_tmp = sim;

duplicate_pairs = [];

continue_to_check = true;

while continue_to_check

    continue_to_check = false;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check for coincident spike fraction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [val_max, arg_max] = max(coin_frac_tmp(:));
    if val_max >= coin_thr
        %I corresponds to the component to which coin_frac
        %is normalized for the maximum value
        [I,J] = ind2sub(size(coin_frac_tmp),arg_max);
    else
        I = [];J = [];
        continue_to_check = false;
    end
    
    if ~isempty(I)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check for STA similarity
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if sim_tmp(I,J) >= sim_thr
            %A duplicate pair is found
            duplicate_pairs = [duplicate_pairs;I J];
            %set corresponding rows and cols in criteria matrices to zero,
            %in order to avoid the participation of the duplicate in
            %subsequent checks:
            sim_tmp(I,:) = 0;
            sim_tmp(:,I) = 0;
            coin_frac_tmp(I,:) = 0;
            coin_frac_tmp(:,I) = 0;
            continue_to_check = true;
        else
            %this was not a duplicate, but we do not want to touch this
            %pairwise combination again, hence we set the respective
            %criteria values to zero:
            sim_tmp(I,J) = 0;
            sim_tmp(J,I) = 0;
            coin_frac_tmp(I,J) = 0;
            coin_frac_tmp(J,I) = 0;
            continue_to_check = true;
        end
    
    end    
end

end