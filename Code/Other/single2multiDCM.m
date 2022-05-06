function [copypath] = single2multiDCM(bpcase)
%single2multiDCM Some MYRIAD cases have dicoms in single file format.
%Rewrite to multiple file format. ** Doesn't seem to work with dcm2ana
%function to convert it to nifti, maybe missing some metadata?
%
% - bpcase: folder containing dicom
%
% Author: Yong Ik Jeong
% Date: 2021-11-03

dcmdir = dir(fullfile(bpcase, '*.dcm'));
dcmpath = fullfile(bpcase, dcmdir(1).name);

img = dicomread(dcmpath);
hdr = dicominfo(dcmpath);

bpcasetmp = [bpcase '_temp'];
movefile(bpcase, bpcasetmp);
mkdir(bpcase);

[nr,nc,ns,nt] = size(img);

count = 0;
for s = 1:ns
    for t = 1:nt
        count = count + 1;
        tmphdr = hdr;
        tmphdr.PerFrameFunctionalGroupsSequence = [];
        tmphdr.PerFrameFunctionalGroupsSequence.(['Item_' num2str(count)]) = hdr.PerFrameFunctionalGroupsSequence.(['Item_' num2str(count)]);
        
        perframehdr = hdr.PerFrameFunctionalGroupsSequence.(['Item_' num2str(count)]).Private_2005_140f.Item_1;
        % PerFrameFunctionalGroupsSequence.Item_1.Private_2005_140f.Item_1.ImagePositionPatient
        % is different than
        % PerFrameFunctionalGroupsSequence.Item_1.PlanePositionSequence.Item_1.ImagePositionPatient
        % why??????
        perframefields = fieldnames(perframehdr);
        perframevalues = struct2cell(perframehdr);
        
        for i = 1:length(perframefields)
            tmphdr.(perframefields{i}) = perframevalues{i};
        end
        
        tmphdr.RepetitionTime = hdr.SharedFunctionalGroupsSequence.Item_1.MRTimingAndRelatedParametersSequence.Item_1.RepetitionTime;
        
        imgpos = hdr.PerFrameFunctionalGroupsSequence.(['Item_' num2str(count)]).PlanePositionSequence.Item_1.ImagePositionPatient;
        imgorient = hdr.PerFrameFunctionalGroupsSequence.(['Item_' num2str(count)]).PlaneOrientationSequence.Item_1.ImageOrientationPatient;
        tmphdr.SliceLocation = imgpos(3);
        tmphdr.ImagePositionPatient = imgpos;
        tmphdr.ImageOrientationPatient = imgorient;
        
        dicomwrite(squeeze(img(:,:,s,t)), fullfile(bpcase, [dcmdir(1).name(1:end-4) '_' num2str(count) '.dcm']), tmphdr, 'CreateMode', 'copy');
    end
end
        