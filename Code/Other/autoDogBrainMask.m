function brainmask = autoDogBrainMask(t1pre)

headmask = t1pre ~= 0;
headmask = imfill(headmask,'holes');
cc = bwconncomp(headmask);
pixellist = cc.PixelIdxList;
sizelist = cellfun(@length,pixellist);
[~,maxind] = max(sizelist);
maxcluster = pixellist{maxind};
headmask = zeros(size(headmask));
headmask(maxcluster) = 1;

stat = regionprops(headmask,'centroid');
ccol = round(stat.Centroid(1));
crow = round(stat.Centroid(2));

tempmask = t1pre >= 0 & t1pre < 400;
tempmask = tempmask.*headmask;
tempmask = imfill(tempmask,'holes');

ctop = find(headmask(:,ccol),1,'first');
cbot = find(headmask(:,ccol),1,'last');

[meshc, meshr] = meshgrid(1:size(tempmask,1),1:size(tempmask,2));
maxarea = 0;
for row = ctop:cbot
    if ~tempmask(row,ccol)
        grow = 1;
        ratio = 8/7;
        r2 = 3;
        r1 = r2*ratio;
        while grow
            ellipsemask = ((meshc-ccol)/r1).^2 + ((meshr-row)/r2).^2 <= 1;
            if isempty(find(tempmask.*ellipsemask,1))
                r2 = r2 + 1;
                r1 = r2*ratio;
            else
                grow = 0;
                area = bwarea(ellipsemask);
            end
        end
        if area > maxarea
            maxarea = area;
            brainmask = ellipsemask;
            brainloc = [row ccol];
        end
    end
end

end
