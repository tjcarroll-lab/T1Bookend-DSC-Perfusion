function icad_uc_nii2dcm(bpnii,bpdcm,bphdr)
%nii2dcm Convert Nifti files to dicom after motion correction
%   Reads in .nii files from bpnii and writes out to dicoms in bpdcm using
%   header info from dicoms in bphdr
%
%   bpnii: folder containing .nii files
%       ex. '\MYRIAD MOTION2\01-006\perf\mc'
%   bpdcm: folder to write dicoms to
%       ex. '\MYRIAD MOTION2\01-006\perf\mc\dcm'
%   bphdr: folder with original dicoms to get header info
%       ex. '\MYRIAD TO PROCESS\01-006\t01\P001\ep2d_perf'
%
% Author: Yong Ik Jeong
% Date: 2021-11-09

% bpparts = split(bpnii, '\');
% caseid = bpparts{end-2};

% Read in Nifti files
niifiles = ls(fullfile(bpnii, '*.nii'));
niifiles = fullfile(bpnii, cellstr(niifiles));
P = spm_vol(char(niifiles));
data = spm_read_vols(P);

% Data from nifti files are in a different orientation so rotate it
data = imrotate(data, 90);
data = permute(data, [1 2 4 3]);
%data = uint16(data);

% Get header info and write to dicom
hdrdir = dir([bphdr '\*.dcm']);
dcmnames = {hdrdir.name};
dcmnames = sortFiles(dcmnames, '-_.');

hdrs = cell(length(hdrdir), 1);
for i = 1:length(hdrdir)
    hdrs{i} = dicominfo(fullfile(bphdr, dcmnames{i}));
end

count = 0;
for s = 1:size(data, 4) % By slice
    for t = 1:size(data, 3) % By time
        count = count + 1;
        tmpdata = uint16((data(:,:,t,s)-P(t).pinfo(2)) / P(t).pinfo(1));
        dicomwrite(tmpdata, [bpdcm '\' dcmnames{count}], hdrs{count}, 'CreateMode', 'copy', 'WritePrivate', true);
    end
end

end

