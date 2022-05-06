function [ perfstruct ] = getMicrospheres_day2( filepath )
%getMicrospheres Get microsphere cbf values from excel sheet
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2017-07-30

[~,~,raw] = xlsread(filepath,'raw');
count = 0;

for ii = 2:1:size(raw,2)
    casenum = raw{1,ii};
    numbers = raw(2:end,1);
    cbfpre = raw(2:end,ii);
    cond1 = cellfun(@isnan,cbfpre);
    cond2 = cellfun(@isempty,cbfpre);
    cond3 = cell2mat(cellfun(@(x) x==0,cbfpre,'UniformOutput',false));
    numbers(cond1 | cond2 | cond3) = [];
    cbfpre(cond1 | cond2 | cond3) = [];
    
    count = count + 1;
    perfstruct(count).casenum = casenum;
    perfstruct(count).cbfpre = [numbers, cbfpre];
    
end
    

end

