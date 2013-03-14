function [STA] = GetSTA(data,time,sr,show)
%GetSTA
%
% Input:
%
% data: (N_ROW x N_COL x N_SAMPLES) array
% time: list of time stamps in ms
% sr: Sampling rate in kHz

[N_ROW,N_COL,N_SAMPLES] = size(data);


spk_tr = round(time * sr);
pre = round(0.5 * sr);
post = round(0.5 * sr);
f_tot = pre + post + 1;


%take only those spikes for which the desired window is entirely contained
%in the data:
spk_tr = spk_tr( ((spk_tr - pre) >= 1) & ((spk_tr + post) <= N_SAMPLES) );

N_spks = length(spk_tr);

%X = reshape(data, [N_ROW*N_COL N_SAMPLES]);

spks = zeros(N_spks,N_ROW,N_COL,f_tot);

for i = 1:N_spks
    spks(i,:,:,:) = data(:,:,spk_tr(i)-pre:spk_tr(i)+post);
end

STA = squeeze(mean(spks,1));
%Get position of maximum:
[row_max,col_max] = find(max(abs(STA),[],3) == max(max(max(abs(STA)))));

if show
    figure;
    for r = 1:3
        for c = 1:3
            subplot(3,3,sub2ind([3 3],c,r));
            if ((row_max-1+r)>0) && ((col_max-1+c)>0) && ...
                    ((row_max-1+r)<=N_ROW) && ((col_max-1+c)<=N_COL)
                plot(squeeze(spks(:,row_max-1+r,col_max-1+c,:))','Color',[0.75 0.75 0.75]);
                hold on;plot(squeeze(STA(row_max-1+r,col_max-1+c,:)),'g')
            end
        end
    end
    figure;
    for i = 1:size(STA,3);
        subplot(ceil(sqrt(size(STA,3))),ceil(sqrt(size(STA,3))),i);
        imagesc(STA(:,:,i));
    end
end





end