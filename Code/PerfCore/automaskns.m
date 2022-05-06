function mask = automaskns(images,threshold)

%%========================================================================%
%%                    making mask without scalp 
%%========================================================================%
% Input  1) images (2D or 3D-matrix) : stack of images
%        2) threshold (scalar) : zero less than threshold*mean 
%           if threshold doesn't exist, default value is 0.7
%%========================================================================%
% Output 1) mask (2D-matrix) : voxel in ROI = 1 otherwise 0
%%========================================================================%
%                                                                02-22-2004 

% by default, threshold is 0.7
if nargin < 2
    threshold = 0.7;    % default value
end

means = sum(sum(images(:,:,1)))/(size(images,1)*size(images,2));
mask = zeros(size(images,1),size(images,2));
cutoffs = threshold*means; % set threshold
mask(find(images(:,:,1) > cutoffs))=1;

