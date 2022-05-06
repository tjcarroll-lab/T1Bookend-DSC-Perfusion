%% This script runs through the motion correction process for ICAD_UC data
% 1. Convert to nifti
% 2. Realign to first time point
% 3. Convert back to dicom
%
% Author: Yong Ik Jeong
% Date: 2021.12.13

%addpath('D:\Users\CarrollLab\Desktop\Yongs Research\Global\MATLAB FILES\MatlabLibrary\spm8');
addpath('D:\Users\CarrollLab\Desktop\Yongs Research\Global\MATLAB FILES\spm12');
addpath('D:\Users\CarrollLab\Desktop\Yongs Research\VascTerr\Code\VT Core');

%% Motion correct perfusion
ptdir = 'D:\Users\CarrollLab\Desktop\Yongs Research\ICAD_UC\ICAD_UC_001_12072021\PERFUSION\P001\';
perfdir = fullfile(ptdir, 'ep2d_perf');
t1predir = fullfile(ptdir, 'LL_EPI_PRE');
t1postdir = fullfile(ptdir, 'LL_EPI_POST');
perfniidir = fullfile(ptdir, 'niifiles', 'perf');
t1preniidir = fullfile(ptdir, 'niifiles', 't1pre');
t1postniidir = fullfile(ptdir, 'niifiles', 't1post');

tmplist = {perfniidir, t1preniidir, t1postniidir};
for i = 1:numel(tmplist)
    if ~exist(tmplist{i}, 'dir')
        mkdir(tmplist{i});
    end
end

if ~exist(fullfile(perfniidir, 'perf001.nii'), 'file')
    dcm2nii(perfdir, fullfile(perfniidir,'perf'));
end

dcm2nii(t1predir, fullfile(t1preniidir, 't1pre'));
dcm2nii(t1postdir, fullfile(t1postniidir, 't1post'));
fprintf('Convert to nii done\n');

% Something's wrong with the T1 scans (both pre, post), so SPM runs into an
% error when it tries to realign
%outdir = mcMyriad(t1preniidir);
%outdir = mcMyriad(t1postniidir);

mcdir = mcMyriad(perfniidir);
fprintf('Realignment done\n');

bpnii = mcdir;
bpdcm = fullfile(mcdir, 'dcm');
if ~exist(bpdcm, 'dir')
    mkdir(bpdcm);
end
bphdr = perfdir;
icad_uc_nii2dcm(bpnii, bpdcm, bphdr);
movefile(bphdr, fullfile(bphdr, '..', 'orig_ep2d_perf'));   
movefile(bpdcm, bphdr);
fprintf('Convert to dicom done\n');

%% Motion correct DCE
leakdir = 'D:\Users\CarrollLab\Desktop\Yongs Research\ICAD_UC\ICAD_UC_001_12072021\LEAKAGE\';
dcedir = fullfile(leakdir, 'ScannerRecon');
dceniidir = fullfile(leakdir, 'niifiles', 'dce');
if ~exist(dceniidir, 'dir')
    mkdir(dceniidir);
end

% Rename files to end with .dcm to work with existing code
dcedirlist = dir(fullfile(dcedir, 'IM*'));
if ~endsWith(dcedirlist(1).name, '.dcm', 'IgnoreCase', true)
    for i = 1:length(dcedirlist)
        dcefilename = fullfile(dcedirlist(i).folder, dcedirlist(i).name);
        movefile(dcefilename, [dcefilename '.dcm']);
    end
end

dcm2nii(dcedir, fullfile(dceniidir, 'dce'));
fprintf('Convert to nii done\n');

mcdir = mcMyriad(dceniidir);
fprintf('Realignment done\n');

bpnii = mcdir;
bpdcm = fullfile(mcdir, 'dcm');
if ~exist(bpdcm, 'dir')
    mkdir(bpdcm);
end
bphdr = dcedir;
icad_uc_nii2dcm(bpnii, bpdcm, bphdr);
movefile(bphdr, fullfile(bphdr, '..', 'orig_ScannerRecon'));   
movefile(bpdcm, bphdr);
fprintf('Convert to dicom done\n');
