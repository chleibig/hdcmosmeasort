function [Xbar] = spatiotemporalembedding(X, Tf)
% [Xbar] = spatiotemporalembedding(X, Tf)
%
% Input
% =====
% 
% X  -    (N,T)-array
% Tf -    Number of samples to embed along the second dimension of X
%
% Output
% ======
% 
% Xbar -  (N*Tf, T)-array, if necessary padded with zeros at the end
%
% Christian Leibig, 16.02.15
%

assert(length(size(X)) == 2, 'X should be a matrix.');

[N, T] = size(X);
%Xbar = zeros(N*Tf,T - mod(T,Tf));
%Xbar = zeros(N*Tf,T);
for t = 1:T-Tf
    Xbar(:,t) = reshape(X(:,t:t+Tf-1)',[N*Tf 1]);
end


end