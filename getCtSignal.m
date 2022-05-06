function [ signal ] = getCtSignal( dsc_slice,mask )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2017-08-09

tp = size(dsc_slice,3);
signal = zeros(1,tp);
for ii = 1:size(dsc_slice,1)
    for jj = 1:size(dsc_slice,2)
        if mask(ii,jj)
            signal = signal + reshape(dsc_slice(ii,jj,:),[1 tp]);
        end
    end
end
signal = signal/sum(mask(:));


end

