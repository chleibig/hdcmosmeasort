function [ S ] = setfieldifnotpresent(S,field,V)
%[ S ] = setfieldifnotpresent(S,field,V) sets the contents of the 
%specified field to the value V, but only if the field is not already
%present
%
% Input
% =====
%
% S         -       1-by-1 struct
% field     -       string
% V         -       anything
%
% christian.leibig@g-node.org, 2015-04-03
%

if isfield(S,field)
else
   S.(field) = V;
end
 
end

