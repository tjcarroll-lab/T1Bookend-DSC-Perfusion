function [outtable] = getROIvalues(img,roi)
%getROIvalues Get mean value of ROI in table format
%   img: image matrix
%   roi: roi matrix
%
% Author: Yong Ik Jeong
% Date: 2020-07-25
% Changelog:
%   - 20200725 YIJ: Initial version


roitable = readtable('..\excel files\ROI Label Collat 31 and after.xlsx');

roiindex = unique(roitable.Index);
roiindex(roiindex == 0) = [];
roivalue = zeros(length(roiindex),1);
roilabel = cell(length(roiindex),1);
for ii = 1:length(roiindex)
    roivalue(ii) = mean(img(roi == roiindex(ii)));
    roilabel{ii} = roitable.Label{roitable.Index == roiindex(ii)};
end
outtable = table(roilabel,roivalue,'VariableNames',{'Label','Value'});

end

