function [SM, Tau] = channel_crosstalk(X, maxlags, SD, d_max)
%Estimate pairwise crosstalk between channels. Currently based
%on correlation coefficient, but consider crosstalk metric used
%by Mads Dyrholm
%crosstalk is only estimated between components i and j if spatial distance
%SD(i,j) between them is maximally d_max (in \mum)

fprintf('Estimating crosstalk...');
[D,T] = size(X);
C_max = zeros(D,D);
Tau_max = zeros(D,D);
t1 = clock;
for i=1:D
    for j=i+1:D
        if SD(i,j) > d_max; C_max(i,j) = 0;
        else
            [c_max,tau_max] = max(abs(xcorr(X(i,:),X(j,:),maxlags,'coeff')));
            C_max(i,j) = c_max;
            Tau_max(i,j) = tau_max;
        end
    end
end
t2 = clock;
fprintf('done in %g sec.\n',etime(t2,t1));

SM = C_max;
Tau = Tau_max;
end