function [  ] = viewDSCsignalv2( dsc_slice, tp, coord )
%viewDSCsignalv2 Summary of this function goes here
%   dsc_slice is a time series of one slice
%
% Author: Yong Ik Jeong
% Date: 2017-12-06

if nargin < 3
    if nargin < 2
        tp = 1;
    end
    figure;
    imshow(dsc_slice(:,:,tp),[],'border','tight');
    [col,row] = ginput();
    row = round(row);
    col = round(col);
%     roiind = find(roi > 0);
%     [row,col] = ind2sub([size(dsc_slice,1) size(dsc_slice,2)],roiind);
else
    row = coord(1);
    col = coord(2);
end

npr = 3;
npc = 3;
[gridc,gridr] = meshgrid([-1 0 1]);
gridc = gridc';
gridr = gridr';

for ii = 1:length(row)
    figure;
    dscary = zeros(npr*npc,size(dsc_slice,3));
    for jj = 1:npr*npc        
        subplot(npr,npc,jj);
        dscary(jj,:) = reshape(dsc_slice(gridr(jj)+row(ii),gridc(jj)+col(ii),:),[1 size(dsc_slice,3)]);
        plot(dscary(jj,:));
        title([num2str(row(ii)+gridr(jj)) ',' num2str(col(ii)+gridc(jj))]);
    end
end


end

