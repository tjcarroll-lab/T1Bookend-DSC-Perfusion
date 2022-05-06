function outdir = mcMyriad(bp)
%mcMyriad Register perfusion series to correct for motion
%   Function uses SPM to register. SPM registers to first volume.
%   Registered files are written to bpcase\mc folder
%
%   bp: folder containing .nii files
%       ex. 'D:\Users\CarrollLab\Downloads\MYRIAD MOTION2\01-006\perf'
%
% Author: Yong Ik Jeong
% Date: 2021-11-09

niifiles = ls(fullfile(bp, '*.nii'));
niifiles = fullfile(bp, cellstr(niifiles));

P = spm_vol(char(niifiles));

clearvars flags1 flags2;

% Realign (calculate transformation matrix)
flags1.quality = 1; % Higher means more precise
flags1.sep = 2; % Sampling distance (2 mm)
flags1.rtm = 0; % Register to mean of images (false; register to first img)
P = spm_realign(P, flags1);

% Reslice (actually resample volumes)
flags2.mask = 1; % Mask out values outside of range of original image
flags2.mean = 0; % Don't write out mean image
flags2.which = 2; % Reslice all volumes
flags2.wrap = [1 1 0]; % Wrapping in x, y direction
spm_reslice(P, flags2);

[filepath,name,ext] = fileparts(niifiles);
rname = strcat('r', name);
rniifiles = strcat(filepath, '\', rname, ext);
if ~exist([bp '\mc'], 'dir')
    mkdir([bp '\mc']);
end
for i = 1:length(rniifiles)
    movefile(rniifiles{i}, [bp '\mc\' rname{i} ext{i}]);
end

%outdir = strcat(bp, '\mc\', rname, ext);
outdir = strcat(bp, '\mc');

end

