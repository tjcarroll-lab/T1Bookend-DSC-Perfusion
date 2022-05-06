function ROIsMaskWM = fManualMaskWM_SS(T1map_pre,T1map_post)

% Get manual mask for white matter from T1 map
% function : ROIsMaskWM = fAutoMaskWM(T1map,Tesla,mask)
%
% Authoer : Jessy Mouannes
%
% description: WM mask is acquired by manual drawing of ROI
%             
% status: stable
%
% versions
%   11-25-2008 (WYS): Initialization


f1 = figure();
imshow(T1map_pre,[0 1200], 'Border', 'tight'); colormap jet; colorbar;
title('T1 map pre')
truesize(f1, [400, 250]);
f2 = figure();
 imshow(T1map_post,[0 1200], 'Border', 'tight'); colormap jet; colorbar;
 title('T1 map post')
 truesize(f2, [400, 250]);
f3 = figure();
imshow((T1map_pre-T1map_post),[-50 100], 'Border', 'tight'); colormap jet; colorbar;
title('T1 Pre - Post')
truesize(f3, [400, 250]);
f4 = figure();
imshow((T1map_post-T1map_pre),[-50 100 ], 'Border', 'tight'); colormap jet; colorbar;
title('T1 Post - Pre')
truesize(f4, [400, 250]);
f5 = figure();
%imshow(T1map_pre,[0 1200]); colormap jet; colorbar;
imshow((T1map_pre-T1map_post),[-50 50], 'Border', 'tight'); colormap jet; colorbar;
title('T1 map pre-post. Draw WM mask')
truesize(f5, [600, 350]);
f6 = figure();
imshow(T1map_pre,[0 1200],'border','tight'); colorbar; colormap jet;
title('T1 map pre');
truesize(f6,[600,350]);
ROIsMaskWM = roipoly;

