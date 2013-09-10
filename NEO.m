function [ Phi ] = NEO( X )
%NEO applies the non-linear-energy operator
%to the signal X and returns it in Phi
%
%Ref.: Kaiser: "On a simple algorithm to calculate the energy of a signal",
%              Int. Conf. on Acoustics, Speech and Signal Processing, 1990
%
% Input
% =====
%
% X - (N x T) array
%
%
% Output
% ======
%
% Phi - (N x T) array
%
% 
% christian.leibig@g-node.org, 21.08.13
%
%


[N,T] = size(X);
Phi = zeros(N,T);

for i = 1:N
    Phi(i,2:end-1) = X(i,2:end-1).^2 - X(i,3:end).*X(i,1:end-2);
end

end

