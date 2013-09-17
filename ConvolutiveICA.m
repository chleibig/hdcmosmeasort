function [S_cica, A_tau, S_noise, A_noise] = ConvolutiveICA(X,L,A,...
                                          sr,d_row,d_col,N_row,N_col,...
                                          d_max,frames_ROI,do_cICA,varargin)
%Documentation goes here ..................
%
% X: Array of dimension (D,T) containing D channels.
% L: order of convolutive mixing model.
% A: (N,D,L+1) mixing matrix identified
%    by preceding dimensionality reduction from N sensors
%    to D channels (e.g. by instantaneous ICA);
% sr: sampling rate in kHz



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotting = 0;
min_skewness = 0.2;
min_corr = 0.15;
M = 0;
max_cluster_size = 2;
max_iter = 5;
maxlags = ceil(sr);
min_no_peaks = 2;
t_s = 0.5;
t_jitter = 1;
coin_thr = 0.5;
thrFactor = 5;

do_cICA

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
        case 'maxlags'
            maxlags = (varargin{i+1});
        case 'plotting'
            plotting = (varargin{i+1});
        case  'min_skewness'
            min_skewness = (varargin{i+1});
        case 'min_corr'
            min_corr = (varargin{i+1});
        case 'M'
            M = (varargin{i+1});
        case 'approach'
            approach = (varargin{i+1});
            if ~(strcmp(approach,'pair') || strcmp(approach,'cluster'))
                error(['Invalid value for approach: ''' approach '''']);
            end
        case 'max_cluster_size'
            max_cluster_size = (varargin{i+1});
        case 'max_iter'
            max_iter = (varargin{i+1});
        case 'sr'
            sr = (varargin{i+1});
        case 'thrFactor'
            thrFactor = (varargin{i+1});
        case 'min_no_peaks'
            min_no_peaks = (varargin{i+1});
        case 't_s'
            t_s = (varargin{i+1});
        case 't_jitter'
            t_jitter = (varargin{i+1});
        case 'coin_thr'
            coin_thr = (varargin{i+1});
        otherwise
            % Hmmm, something wrong with the parameter string
            error(['Unrecognized parameter: ''' varargin{i} '''']);
    end;
  end;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% S_cica sources;
% A_tau (N_sensors,D_sources,L+1) updated mixing matrices

if (size(A,3) - 1) ~= L
    error('Mismatch of order L and number of mixing matrices A.')
end


iteration_no = 0;
touched = struct('IDs',{});
S_noise = [];A_noise = [];
while iteration_no < max_iter
    iteration_no = iteration_no + 1;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Remove components    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1. get indices of components to keep
    
    [keep] = checkfornoisycomponents(X,min_skewness,thrFactor, min_no_peaks,sr,plotting);

    % 2. remove components and store them away:
    S_noise = [S_noise;X(~keep,:)];
    A_noise = cat(2,A_noise,A(:,~keep,:));
    X = X(keep,:);
    A = A(:,keep,:);
    
    % 3. Adapt touched combination indices
    
    %touched is a struct array with touched(i).IDs containing those ID 
    %combinations that were already processed with convolutive ICA and need
    %not be touched again. When e.g. a cICA component is skipped (due to
    %noise) it makes sense to touch the remaining combination j again because
    %it might have absorbed signal energy from the skipped component(s).
    %hence the respective entry touched(j).IDs is removed.
    
    %skipped IDs:
    delIDs = find(~keep);   
    %keep only combinations for which no member was deleted:
    touched = touched(cellfun(@(x) ~any(ismember(x,delIDs)),{touched.IDs}));
    %adapt remaining indices:
    for i=1:length(touched)
        for j=1:length(touched(i).IDs)
            oldID = touched(i).IDs(j);
            touched(i).IDs(j) = oldID - nnz(delIDs < oldID);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimate crosstalk
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    SD = SpatialDistance(A, d_row, d_col, N_row, N_col);

    SM = channel_crosstalk_wiener(X(:,frames_ROI), maxlags, SD, d_max);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Group channels based on crosstalk
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [T] = hierarchical_clustering(SM,min_corr,'plotting',plotting);
    %counts(i): total # of members belonging to cluster i:
    counts = arrayfun(@(x)(sum(T == x)),1:max(T));
    fprintf('Identified %g clusters showing crosstalk.\n',...
        length(nonzeros(counts - 1)));
    fprintf('The largest cluster contains %g channels.\n',max(counts));
    %care only about clusters that have more than one member:
    cluster_ids = find(counts>=2);
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reduce cluster size
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Idea: if maximum cluster size (due to computational load) is exceeded
    %unmix subcluster only, this is independent from all other clusters, 
    %within each bigger cluster, we can only unmix the remaining components
    %in a sequential manner (not in parallel); idea: increase
    %minimum similarity until maximum sub cluster size <= max_cluster_size
    if max(counts) > max_cluster_size
        fprintf('Only the %g channels showing the strongest\n',max_cluster_size);
        fprintf('crosstalk within each cluster will be unmixed.\n');

        for i = 1:length(cluster_ids)
            cl_i = find(T == cluster_ids(i))';
            if (length(cl_i) > max_cluster_size)
                d_sm = 0;
                sm_cl_i = SM(cl_i,cl_i);
                cl_i_tmp = cl_i;
                while length(cl_i_tmp) > max_cluster_size
                    d_sm = d_sm + 0.01;
                    t = hierarchical_clustering(sm_cl_i,...
                        min_corr + d_sm,'plotting',0);
                    t_counts = arrayfun(@(x)(sum(t == x)),1:max(t));
                    %take the biggest remaining subcluster:
                    [counts_max,ind_max] = max(t_counts);
                    cl_i_tmp = cl_i(t == ind_max);
                    if length(cl_i_tmp) <= max_cluster_size
                        cl_i = cl_i_tmp;
                    end
                end
                clear cl_i_tmp;
                %cl_i is the subcluster to unmix, put each of the remaining
                %components in a separate cluster (because due to single
                %linkage clustering they do not necessarily have something in
                %common):
                cl_rem = setdiff(find(T == cluster_ids(i)),cl_i);
                for k = 1:length(cl_rem)
                    T(cl_rem(k)) = max(T) + 1;
                end
                clear cl_rem
            end
        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Skip clusters that were already touched
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %all clusters have to be compared against all touched combinations
    %for equality:
    equality = zeros(length(cluster_ids),length(touched));
    for i = 1:length(cluster_ids)
        for j = 1:length(touched)
            c = find(T == cluster_ids(i));
            equality(i,j) = isequal(sort(c(:)),sort(touched(j).IDs(:)));
        end
    end
    N_total = length(cluster_ids);
    cluster_ids_touched = cluster_ids(any(equality'));
    N_touched = length(cluster_ids_touched);
    %Put each component of the already touched ones in a separate cluster
    %to guarantee, that T only contains clusters with more than one member
    %that we actually want to provide as input to cICA:
    for i = 1:length(cluster_ids_touched)
        cl_i_touched = find(T == cluster_ids_touched(i));
        for k = 1:length(cl_i_touched)
            %if k == 1: the first component is left as only member with
            %ID T(cl_i_touched(1))
            if k > 1 
                T(cl_i_touched(k)) = max(T) + 1;
            end
        end
        clear cl_i_touched
    end

    %cluster_ids = cluster_ids(~any(equality'));
    
    fprintf('%g of the %g total clusters were already touched\n',...
        N_touched,N_total); 
    
    % Keep only clusters with more than one member:
    counts = arrayfun(@(x)(sum(T == x)),1:max(T));
    cluster_ids = find(counts>=2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Crosstalk before cICA step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if plotting
        %Crosstalk
        figure;title('Clusterwise crosstalk');
        pltsize = ceil(sqrt(max(T)));
        for i = 1:max(T)
            subplot(pltsize,pltsize,i);plot(X(find(T == i),:)');
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Convolutive ICA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    if ~isempty(cluster_ids) && iteration_no < max_iter && do_cICA
        fprintf('Iteration %g...\n',iteration_no);
        switch approach
            case 'cluster'
                for i=1:length(cluster_ids)
                    fprintf('CICAAR for unmixing cluster %g,',i);
                    fprintf(' containing channels:\n');
                    cl_i = find(T == cluster_ids(i))'
                    %this is the most time consuming step:
                    tic;
                    [e,ll,bic,invA0,Atau,Hlambda] = cicaarpro(X(cl_i,frames_ROI),L,M);
                    toc;
                    fprintf('BIC=%f, e=%f\n',bic,e);
                    %Mixing matrix update:
                    A_new = cat(3,inv(invA0),Atau);
                    A_old = A(:,cl_i,:);
                    A(:,cl_i,:) = 0;
                    %the following could be vectorized:
                    for tau = 1:L+1
                        for k = 1:tau
                            j = tau + 1 - k;
                            A(:,cl_i,tau) = A(:,cl_i,tau) + ...
                                            A_old(:,:,k)*A_new(:,:,j);
                        end
                    end
                    X(cl_i,:) = cicaarprosep(invA0,Atau,X(cl_i,:));
                    %X = cicaarprosep(pinv(A(:,:,1)),A(:,:,2:end),X_org); 
                    %append current cluster to list of touched component
                    %combinations:
                    touched(length(touched)+1).IDs = cl_i;
                end
                
            case 'pair' %maybe conceptually wrong and should be removed
                [I,J] = find(abs(SM) >= min_corr);
                M_pairs = length(I);
                invA0_tmp = eye(size(X,1),size(X,1));
                Atau_tmp = zeros(size(X,1),size(X,1),L);
                for i=1:M_pairs
                    fprintf('CICAAR for unmixing pair %g ...',i);
                    %this is the most time consuming step:
                    tic;
                    [e,ll,bic,invA0,Atau,Hlambda] = cicaarpro(X([I(i) J(i)],:),L,M);
                    toc;
                    invA0_tmp([I(i) J(i)],[I(i) J(i)]) = invA0;
                    Atau_tmp([I(i) J(i)],[I(i) J(i)],:) = Atau;
                    fprintf('BIC=%f, e=%f\n',bic,e);
                    A_new = cat(3,inv(invA0),Atau);
                    A_old = A(:,[I(i) J(i)],:);
                    A(:,[I(i) J(i)],:) = 0;
                    for tau = 1:L+1
                        for k = 1:tau
                            j = tau + 1 - k;
                            A(:,[I(i) J(i)],tau) = A(:,[I(i) J(i)],tau)+...
                                                 A_old(:,:,k)*A_new(:,:,j);
                        end
                    end
                    touched(length(touched)+1).IDs = [I(i) J(i)];
                end
                %X = cicaarprosep(pinv(A(:,:,1)),A(:,:,2:end),X_org);
                %the following unmixing was observed to lead to
                %instabilities:
                X = cicaarprosep(invA0_tmp,Atau_tmp,X);
        end

        if plotting
            %Crosstalk
            figure;title('Remaining Clusterwise crosstalk');
            pltsize = ceil(sqrt(max(T)));
            for i = 1:max(T)
                subplot(pltsize,pltsize,i);plot(X(find(T == i),:)');
            end
        end
        
    else
        S_cica = X;
        A_tau = A;
        if iteration_no < max_iter;
            fprintf('Converged in %g iterations.\n',iteration_no -1);
        else
            fprintf('Maximum number of iterations reached (%g)\n',max_iter);
        end
        
        return;
    end
    

end
 

end



