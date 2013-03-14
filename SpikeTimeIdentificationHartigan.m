function [units] = SpikeTimeIdentificationHartigan(X,sr,show)
% Perform spike time identification on each component
% of X (dims x samples) individually via 
% determining threshold crossing events. Potentially presence
%of more than one neuron on a single component is assessed
%via testing for multimodality with hardigans dip test and
%in case of multimodality setting the threshold to the LOWEST or HIGHEST??
%local minimum

N_comp = size(X,1);
%correct skewness:
X(skewness(X') > 0,:) = -1 * X(skewness(X') > 0,:);

%For Hardigans diptest:
dip = zeros(N_comp,1);
p = zeros(N_comp,1);
xlow = zeros(N_comp,1);
xup = zeros(N_comp,1);
nboot = 500; %number of bootstrap samples drawn from uniform pdf
sign_lev = 0.05;

if(show);fX=figure;fH=figure;splt_size = ceil(sqrt(N_comp));end
units = struct('time',{},'amplitude',{});
for i=1:N_comp;
    if(show);
        figure(fX);subplot(splt_size,splt_size,i);plot(X(i,:));hold on;
    end
    [indices,pos_amplitudes] = find_peaks(-X(i,:),...
                                    5*median(abs(X(i,:))/0.6745),ceil(sr));
    amplitudes = -1*pos_amplitudes;
    
    if ~isempty(amplitudes)
        xpdf = sort(amplitudes);
        [dip(i), p(i), xlow(i), xup(i)] = HartigansDipSignifTest(xpdf, nboot);
        [cts, bin_ctrs] = hist(xpdf,floor(sqrt(length(xpdf))));
        if show
            %all threshold crossings:
            figure(fX);subplot(splt_size,splt_size,i);
            plot(indices,X(i,indices),'rx','LineStyle','none');            
            %hist with dip test results
            figure(fH);subplot(splt_size,splt_size,i);hold on;
            bar(bin_ctrs,cts);
            title(['dip=',num2str(dip(i),3), ', p=',num2str(p(i),3)]);
        end
        if p(i) < sign_lev %adapt threshold
            [val,ind] = min(cts);
            thr_loc_min = bin_ctrs(ind);
            indices = indices(amplitudes < thr_loc_min);
            units(i).time = indices/sr;
            units(i).amplitude = X(i,indices);
            if show
                figure(fH);plot([thr_loc_min thr_loc_min], [0 max(cts)]);
                figure(fX);plot(indices,X(i,indices),'go','LineStyle','none');             
            end
            
        else
            units(i).time = indices/sr;
            units(i).amplitude = X(i,indices);
        end
    else
        units(i).time = [];
        units(i).amplitude = [];
    
    end
end



end
