function [ dsc1, aifmap ] = getAIFmap( dsc, aifpos )
%getAIFmap Get single timepoint dsc stack and aif position map
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2017-05-24

%aifpos = ROIs.positions.AIF;
%dsc = ROIs.dsc_stack;

nslice = length(dsc);
ntime = size(dsc{1},3);
totalsize = [size(dsc{1},1) size(dsc{1},2) nslice];

aifmap = zeros(totalsize);
for ii = 1:length(aifpos)
    aifmap(aifpos{ii}(1),aifpos{ii}(2),aifpos{ii}(3)) = 1;
end

dsc1 = zeros(totalsize);
for ii = 1:nslice
    dsc1(:,:,ii) = dsc{ii}(:,:,1);
end

figure;
for ii = 1:length(aifpos)
    subplot(round(sqrt(length(aifpos))),ceil(sqrt(length(aifpos))),ii);
    plot(smooth(squeeze(dsc{aifpos{ii}(3)}(aifpos{ii}(1),aifpos{ii}(2),:)),5,'sgolay',3));
end

end

