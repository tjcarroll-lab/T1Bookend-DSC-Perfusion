function [avg sigma] = fast_roi(true_roi, image)
%-------------------------------------------------------------------------%
% Format: [avg sigma] = fast_roi(true_roi, image)
% 
% Input:
% true_roi: The roi(s) mask to use for roi information
% image: the image(s) of which roi is made from
% 
% Output:
% avg: mean value of roi
% sigma: STD value of roi
% 
% This function adapted (to cover image volume) from Dr. Carroll's version 
% takes the mean and std flow value given roi(s) of particular slice(s).
% 
% Author(of revision): Maulin Shah 
% Date: 10/13/2006
%-------------------------------------------------------------------------%

[x,y,z] = size(image);
N_voxels = 0.0; 
total_image = 0.0;
for i = 1:x
    for j = 1:y
        for k = 1:z
            if ((true_roi(i,j,k) > 0.0 ) && (image(i,j,k) > 0.0)) %Jessy (3/18/2008): exclude zero-voxels
                N_voxels = N_voxels + 1.0;
                total_image = total_image + image(i,j,k);
            end
        end
    end
end
avg = total_image/N_voxels;

std_sum = 0.0;
for i = 1:x
    for j = 1:y
        for k = 1:z
            if (( true_roi(i,j,k) > 0.0) && (image(i,j,k) > 0.0)) %Jessy (3/18/2008): exclude zero-voxels
                std_sum = std_sum + (double(image(i,j,k))-avg)^2;
            end
        end
    end
end

sigma = sqrt(std_sum/N_voxels);

return