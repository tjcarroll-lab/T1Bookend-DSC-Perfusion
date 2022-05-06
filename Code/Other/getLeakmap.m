function [ dsc1, enhancemap, sumsignals ] = getLeakmap( dsc )
%getLeakmap 
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2017-05-24

%dsc = ROIs.dsc_stack;

nslice = length(dsc);
ntime = size(dsc{1},3);
totalsize = [size(dsc{1},1) size(dsc{1},2) nslice];

enhancemap = zeros(totalsize);
for ii = 1:nslice
    tmpmat = reshape(dsc{ii},totalsize(1)*totalsize(2),ntime)';
    mean1 = mean(tmpmat(1:50,:));
    mean2 = mean(tmpmat(150:200,:));
    tmpmap = zeros(totalsize(1),totalsize(2));
    tmpmap(mean2 > mean1) = 1;
    enhancemap(:,:,ii) = tmpmap;
end

dsc1 = zeros(totalsize);
brainmask = zeros(totalsize);
f1 = figure;
for ii = 1:nslice
    dsc1(:,:,ii) = dsc{ii}(:,:,1);
    figure(f1);
    imshow(dsc1(:,:,ii),[]);
    brainmask(:,:,ii) = roipoly;
end
close(f1);
enhancemap = enhancemap.*brainmask;

totalpixels = sum(enhancemap(:));
sumsignals = zeros(1,ntime);
for ii = 1:nslice
    tmpdsc = dsc{ii};
    tmpmat = reshape(tmpdsc,totalsize(1)*totalsize(2),ntime);
    tmpmat2 = reshape(enhancemap(:,:,ii),totalsize(1)*totalsize(2),1);
    sumsignals = sumsignals + sum(tmpmat(logical(tmpmat2),:));
end
sumsignals = sumsignals / totalpixels;

end

