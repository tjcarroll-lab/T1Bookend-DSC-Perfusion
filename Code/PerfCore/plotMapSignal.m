function [] = plotMapSignal( dsc, mask )
%plotMapSignal Get signal curves of AIF/Vein map mask 
%   dsc: 4D DSC matrix
%   mask: 3D mask of AIF/Vein locations
%
% Author: Yong Ik Jeong
% Date: 2021-02-03

%aifpos = ROIs.positions.AIF;
%dsc = ROIs.dsc_stack;

nslice = size(dsc,4);
totalsize = [size(dsc,1) size(dsc,2) nslice];

aifpos = find(mask);
[rr,cc,kk] = ind2sub(size(mask),aifpos);

figure;
for ii = 1:length(aifpos)
    subplot(round(sqrt(length(aifpos))),ceil(sqrt(length(aifpos))),ii);
    plot(smooth(double(reshape(dsc(rr(ii), cc(ii), :, kk(ii)), [1 size(dsc,3)])), 5, 'sgolay', 3));
    title(sprintf('r%d c%d s%d', rr(ii), cc(ii), kk(ii)));
end

end