function [ aifmap ] = getAIFmapv2( dsc, aifpos )
%getAIFmap Get single timepoint dsc stack and aif position map
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2017-12-06

%aifpos = ROIs.positions.AIF;
%dsc = ROIs.dsc_stack;

nslice = size(dsc,4);
totalsize = [size(dsc,1) size(dsc,2) nslice];

aifmap = zeros(totalsize);
for ii = 1:length(aifpos)
    aifmap(aifpos{ii}(1),aifpos{ii}(2),aifpos{ii}(3)) = 1;
end

figure;
for ii = 1:length(aifpos)
    subplot(round(sqrt(length(aifpos))),ceil(sqrt(length(aifpos))),ii);
    plot(smooth(double(reshape(dsc(aifpos{ii}(1),aifpos{ii}(2),:,aifpos{ii}(3)),[1 size(dsc,3)])),5,'sgolay',3));
end

end

