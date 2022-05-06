function [ roi_stack ] = draw_mcrsphrROIs( dsc_stack,num_rois_per_slc )
%draw_mcrsphrROIs draw ROIs for microsphere comparison 
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2017-07-05

f1 = figure;
roicount = 0;
%roi_stack = zeros(size(dsc_stack{1},1),size(dsc_stack{1},2),length(dsc_stack));
roi_stack = zeros(size(dsc_stack));
%for ii = 1:length(dsc_stack)
for ii = 1:size(dsc_stack,3)
    %dsc_slc = dsc_stack{ii}(:,:,1);
    dsc_slc = dsc_stack(:,:,ii);
    num_rois = num_rois_per_slc(ii);
    roi_slc = roi_stack(:,:,ii);
    
    if num_rois > 0
        for jj = 1:num_rois
            yv(dsc_slc,[],'overlay',roi_slc);
            tmpf = gcf;
            figure(f1);imshow(dsc_slc,[]);
            tmproi = roipoly;
            roinum = input('ROI num: ');
            roi_slc(tmproi) = roinum;
            close(tmpf);
        end
    end
    
    roi_stack(:,:,ii) = roi_slc;
end

end
        
        
        