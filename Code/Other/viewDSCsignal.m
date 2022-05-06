function [ dscary ] = viewDSCsignal( dsc_slice )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

figure;
imshow(dsc_slice(:,:,1),[]);
roi = roipoly;

roiind = find(roi > 0);
[row,col] = ind2sub([size(dsc_slice,1) size(dsc_slice,2)],roiind);

dscary = zeros(length(roiind),size(dsc_slice,3));
plotcount = 0;
figure;
for ii = 1:length(roiind)
    plotcount = plotcount + 1;
    if plotcount > 9
        figure;
        plotcount = 1;
    end
    subplot(3,3,plotcount);
    dscary(ii,:) = reshape(dsc_slice(row(ii),col(ii),:),[1 size(dsc_slice,3)]);
    plot(dscary(ii,:));
    title([num2str(row(ii)) ',' num2str(col(ii))]);
end


end

