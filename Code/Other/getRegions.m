function [ regions ] = getRegions( filepath )
%getRegions Get region labels and hemispheric numbers
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2017-07-30

[~,~,raw] = xlsread(filepath);
count = 0;
for ii = 1:2:size(raw,2)
    casenum = raw{1,ii};
    regionlist = [];
    hemislist = [];
    for jj = 2:size(raw,1)
        regionname = raw{jj,ii};
        regionind = raw{jj,ii+1};
        if ~sum(isnan(regionname)) && ~sum(isnan(regionind))
        tmpcell = regexpi(regionname,'(\w+)\s+(r|l)\s+(.+)','tokens');
        regionloc = tmpcell{1}{1};
        hemisphere = tmpcell{1}{2};
        sliceloc = tmpcell{1}{3};
        regionlist = [regionlist; {regionname,regionloc,hemisphere,sliceloc,regionind}];
        if isempty(hemislist)
            hemislist = {sliceloc,hemisphere,regionind};
        else
            tmpind = strcmpi(sliceloc,hemislist(:,1)) & strcmpi(hemisphere,hemislist(:,2));
            if sum(tmpind)
                hemislist{tmpind,3} = [hemislist{tmpind,3},regionind];
            else
                hemislist = [hemislist; {sliceloc,hemisphere,regionind}];
            end
        end
        end
    end
    count = count + 1;
    regions(count).casenum = casenum;
    regions(count).regionlist = regionlist;
    regions(count).hemislist = hemislist;
end


end

