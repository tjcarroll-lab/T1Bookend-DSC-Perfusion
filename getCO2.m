function [co2] = getCO2( filepath )
%getCO2 Summary of this function goes here
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2017-11-15

[~,~,raw] = xlsread(filepath);
count = 0;

for ii = 2:size(raw,1)
    casenum = raw{ii,1};
    base = raw{ii,2};
    stress = raw{ii,3};

    count = count + 1;
    co2(count).casenum = casenum;
    co2(count).base = base;
    co2(count).stress = stress;
end



end

