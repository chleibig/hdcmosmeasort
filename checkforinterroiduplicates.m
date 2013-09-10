function [duplicate_pairs, nDupl, mDupl] = checkforinterroiduplicates(nUnits,mUnits, sr, ...
                                           data, t_s, t_jitter, coin_thr,...
                                           sim_thr, show, interactive)
%checkforinterroiduplicates tests all units in nUnits pairwisely against all units
%in mUnits whether they are duplicates or not based on the following criteria:
%
% 1. They have more than coin_thr spikes in common
%
% 2. Their STAs (calculated on data) are more similar than sim_thr
%    (similarity is measured via their normalized scalar product)
%
% INPUT:
% 
% nUnits: array of structs with at least the field "time"
% mUnits: the same as nUnits, but containing DIFFERENT UNITS!
% sr: sampling rate in kHz
% data - (N_ROW x N_COL x N_FRAMES) array 
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
%                  were found to be duplicates; j-th duplicate is made up
%                  by nUnits(duplicate_pairs(j,1)) and
%                  mUnits(duplicate_pairs(j,2))
%
% 
% lessSpikes: binary mask of the same format as duplicate_pairs. The
%             unit that contributes with less spikes is indicated as follows:
%             lessSpikes(j,:) = [true false] <-> nUnits(duplicate_pairs(j,1))
%             lessSpikes(j,:) = [false true] <-> mUnits(duplicate_pairs(j,2))
%             
% nDupl - vector of indices specifying duplicate positions in vector nUnits
% mDupl - vector of indices specifying duplicate positions in vector mUnits

%
%
% author: Christian Leibig, 19.03.13, adapted 22.07.13

%procedure: we check all pairs for the most evident duplicate pair, remove
%           that one and then do the same for the rest and so on until
%           we do not find any more duplicate
%





N = length(nUnits);
M = length(mUnits);

%first compute all the N*M pairwise criteria:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coincident spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coin_frac = zeros(N,M);
for n=1:N
    for m=1:M
        tn = nUnits(n).time;
        tm = mUnits(m).time;
        [no_coin,unused1,unused2] = spiketrainalignment(tn,tm,sr,t_s,t_jitter);
        %choose normalization which maximizes coin_frac
        if length(tn) >= length(tm)
            coin_frac(n,m) = no_coin/length(tm);
        else
            coin_frac(n,m) = no_coin/length(tn);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STA similarities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sim = zeros(N,M);
for n=1:N
    for m=1:M
        %Calculate STAs on the fly:
        STAn = GetSTA(data,nUnits(n).time,sr,0);
        STAm = GetSTA(data,mUnits(m).time,sr,0);
        STAn = reshape(STAn,...
         [size(STAn,1)*size(STAn,2) size(STAn,3)]);
        STAm = reshape(STAm,...
         [size(STAm,1)*size(STAm,2) size(STAm,3)]);
        %temporal alignment on peaks:
        [n_v_p,n_p] = max(max(abs(STAn),[],1));
        [m_v_p, m_p] = max(max(abs(STAm),[],1));
        shift = abs(n_p - m_p);
        if n_p >= m_p
            STAn = STAn(:,1+shift:end);
            STAm = STAm(:,1:end-shift);
        else
            STAm = STAm(:,1+shift:end);
            STAn = STAn(:,1:end-shift);
        end
        STAn = STAn(:);
        STAm = STAm(:);
        %calculate similarity:
        sim(n,m) = STAn'*STAm / sqrt((STAn'*STAn) * (STAm'*STAm));
    end
end


if show
   figure; 
   %this is for debug purposes and to assess suitable threshold parameters
   subplot(2,2,1);imagesc(coin_frac);colorbar;
   title('coincident spikes');
   xlabel('unit index');
   ylabel('unit index');
   %axis square;
   subplot(2,2,3);hist(coin_frac(:));
   xlabel('fraction of coincident spikes');
   ylabel('counts');
   subplot(2,2,2);imagesc(sim);colorbar;
   title('similarities of STAs');
   xlabel('unit index');
   ylabel('unit index');
   %axis square;
   subplot(2,2,4);hist(sim(:));
   xlabel('c_{ij}');
   ylabel('counts');
end

if interactive
   coin_thr = input('Please specify the minimum fraction of coincident spikes for duplicates:');
   sim_thr = input('Please specify the minimum similarity of STAs for duplicates:');
end
    
%*_tmp criteria matrices will be changed in place to keep track of which
%units were already considered as duplicates and should be excluded from
%further comparisons:
coin_frac_tmp = coin_frac; 
sim_tmp = sim;

duplicate_pairs = [];
lessSpikes = logical([]);

continue_to_check = true;

while continue_to_check

    continue_to_check = false;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check for coincident spike fraction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [val_max, arg_max] = max(coin_frac_tmp(:));
    if val_max >= coin_thr
        %I indexes into nUnits, J into mUnits.
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
            %A DUPLICATE IS FOUND
            duplicate_pairs = [duplicate_pairs;I J];
            %check whether nUnits(I) or mUnits(J) contributes with less 
            %spikes and set corresponding entries in criteria matrices to 
            %zero in order to avoid the participation of the duplicate in 
            %subsequent checks:
            if length(nUnits(I).time) >= length(mUnits(J).time)
                %mUnits(J) ought to be removed
                lessSpikes = [lessSpikes; false true];
                %-> Ignore J-th columns: 
                sim_tmp(:,J) = 0;
                coin_frac_tmp(:,J) = 0;
            else
                %nUnits(I) ought to be removed
                lessSpikes = [lessSpikes; true false];
                %-> ignore I-th rows:
                sim_tmp(I,:) = 0;
                coin_frac_tmp(I,:) = 0;
                
            end
            continue_to_check = true;
        else
            %NO DUPLICATE IS FOUND,
            %but we do not want to touch this pairwise combination again,
            %hence we set the respective criteria values to zero:
            sim_tmp(I,J) = 0;
            coin_frac_tmp(I,J) = 0;
            continue_to_check = true;
        end
    
    end    
end

if ~isempty(duplicate_pairs)
    nDupl = duplicate_pairs(lessSpikes(:,1),1);
    mDupl = duplicate_pairs(lessSpikes(:,2),2);
else
    nDupl = [];
    mDupl = [];
end


end