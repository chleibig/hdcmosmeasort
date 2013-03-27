function [ w ] = GetWienerFilterCoeff( input, reference, order )
%WIENERFILTER computes the wiener filter coefficients
% w = [w_0 ... w_order] by minimizing the squared
% error between the filtered input sequence and the
% reference signal (the desired source signal)

% author: Christian Leibig, 26.03.13
%
% implementation follows Wikipedia article on Wiener filter and therein
% "Finite impulse response Wiener filter for discrete series"
%
% Orig. ref.: Wiener, Norbert: "Extrapolation, Interpolation and Smoothing
%             of Stationary Time Series", New York, Wiley, 1949

ac = xcorr(input, order,'unbiased'); %autocorrelation
ac = ac(order+1:2*order+1); %skip negative lags
cc = xcorr(reference,input,order,'unbiased');
cc = cc(order+1:2*order+1)'; %skip negative lags

T = toeplitz(ac);

%solve Wiener-Hopf equations:
w = inv(T) * cc;

end
