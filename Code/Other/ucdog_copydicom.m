function ucdog_copydicom(originpath, newpath)
% Copy and rename dicoms to 1.dcm, 2.dcm, 3.dcm, ... from EPI-LL scans
%
% Author: Yong Ik Jeong
% Date: 2016-03-16

origindir = dir(originpath);
for ii = 3:length(origindir)
    if strfind(origindir(ii).name,'EPI-LL')
        newtemppath = [newpath '\' origindir(ii).name];
        if ~exist(newtemppath,'dir')
            mkdir(newtemppath);
        end
        origindcmpath = [originpath '\' origindir(ii).name];
        origindcmdir = dir([origindcmpath '\*.dcm']);
        for jj = 1:length(origindcmdir)
            copyfile([origindcmpath '\' origindcmdir(jj).name],[newtemppath '\' num2str(jj) '.dcm']);
        end
    end
end
