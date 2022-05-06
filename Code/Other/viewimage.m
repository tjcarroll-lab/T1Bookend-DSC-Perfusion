function viewimage(varargin)
%-------------------------------------------------------------------------%
% Format: viewimage(Image_Vol, scale, maptype)
%   Example: viewimage(qCBF_Vol, [0 120], 'jet');
% 
% Inputs:
% Image_Vol: the volume of images you would like to see
% scale: the upper and lower limit of image scale ex [low high]
% map: the colormap used to view images, default is gray
%
% This function allows you to view a volume of images, with function to
% move back and forth, the left arrow key will move back one slice, while
% right arrow will move forward one slice, exit by pressing space bar
%
% Author: Maulin Shah
% Date: 05/27/2008
%-------------------------------------------------------------------------%
%% Initialize
if nargin >= 1
    Image_Vol = varargin{1};
    if nargin > 1
        scale = varargin{2};
    else
        scale = [];
    end
    if nargin > 2
        maptype = varargin{3};
    else
        maptype = 'gray';
    end
else
    error('Not enough inputs');
end

z = size(Image_Vol,3);
keycode = 0;
slice = 1;

%% Instructions
fprintf('Use forward and back arrows to navigate slices, spacebar to exit:\n');

%% Display Images
while keycode ~= 32

    if keycode == 28
        if slice == 1
            slice = z;
        else
            slice = slice -1;
        end
    end
    if keycode == 29
        if slice == z
            slice = 1;
        else
            slice = slice + 1;
        end
    end

    imshow(Image_Vol(:,:,slice), scale); eval(['colorbar; colormap ' maptype ';']);
    title(['slice #' num2str(slice)]);
    [ignore,ignore2,keycode] = ginput(1);

end

return