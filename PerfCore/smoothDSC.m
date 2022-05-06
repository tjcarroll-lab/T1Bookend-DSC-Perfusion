function [newimages] = smoothDSC(images)
%smoothDSC Return spatially smoothed DSC images
%   xy smoothing and z smoothing
%
% Author: Yong Ik Jeong
% Date: 2020-03-26

newimages = cell(size(images));
for kk = 1:length(images)
    if kk == 1 || kk == length(images)
        newimages{kk} = images{kk};
    else
        newimages{kk} = 0.25*images{kk-1} + 0.5*images{kk} + 0.25*images{kk+1};
    end
    
    mask = automaskns(imfilter(newimages{kk},fspecial('gaussian',[3 3],2)));
    newimages{kk} = fSpatialFilter(newimages{kk},[2 3],mask);
end
