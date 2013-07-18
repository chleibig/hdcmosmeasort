function [ frameIdx ] = time2frame( time, REF_TIME )
%time2frame converts all time [ms] to frame indices based on reference
%           REFERENCE_TIME
%
% 
% Input
% =====
%
% time - a vector of time stamps in ms
% REF_TIME - a vector of reference time stamps in ms
%
% Note: Of course it would be cheaper to simply calculate the frame indices
%       based on the sampling rate. But this requires a constant
%       sampling rate, which might not generally be true.
%
% christian.leibig@g-node.org, 18.07.13
%

frameIdx = zeros(length(time),1);
for i = 1:length(time)
    frameIdx(i) = find( REF_TIME > time(i), 1 ) - 1;
end

end

