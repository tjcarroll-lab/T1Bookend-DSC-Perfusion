function [ sig_out ] = freqfilter( sig_in,type,cutoff )
%freqfilter frequency filter test
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2017-08-23

filter = zeros(size(sig_in));

if strcmpi(type,'low')
    sigma = 100;
    gauss = exp(-.5.*((0:cutoff)/sigma).^2);
    filter(1:cutoff+1) = gauss;
    filter(end-cutoff+1:end) = fliplr(gauss(2:end));
    sig_out = ifft(fft(sig_in).*filter);
else
    sig_out = 0;
end

% sigma = 5;
% gauss = exp(-.5.*((0:(hi-lo))/sigma).^2);
% filter(lo+1:hi+1) = gauss;
% filter(end-hi+1:end-max([0 lo-1])) = fliplr(gauss(2:end));
% sig_out = ifft(fft(sig_in).*filter);

end

