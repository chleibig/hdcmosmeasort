function Z = cicaarwhitensources(Hlambda,S)
% CICAARPROWHITENSOURCES
%
% Synopsis
% ========
%
% Z = cicaarwhitensources(Hlambda,S)
%
% Purpose
% =======
% 
% Whiten sources as estimated by the CICAAR functions CICAARPRO together
% with CICAARPROSEP
% available at http://www.machlea.com/mads/cicaar-pro.html
% 
% the original code is copyright protected as follows:
% 
%  -- Copyright: Mads Dyrholm, all rights reserved --
%     November 2012, Copenhagen, Denmark.
%
%
% The code at hand is authored by Christian Leibig, 2014
%

if isempty(Hlambda)
    Z = S;
    return
end

[M,DH] = size(Hlambda);
[DS,T] = size(S);

assert(DH==DS,['Dimension mismatch between the number of sources and '...
               'whitening filters.']);
           
Z = zeros(DS,T+M);
for d = 1:DS
    for t = 1:T
        Z(d,t+M) = S(d,t) - Z(d,t-1+M:-1:t)*Hlambda(:,d);
    end
end
Z = Z(:,M+1:end);

end