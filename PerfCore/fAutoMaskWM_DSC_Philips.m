function [mask, nslice, index] = fAutoMaskWM_DSC_Philips(path1,path2,maskSS,Nslices,N_meas)
%(T1map,Tesla,mask,dT1)

% Get DSC mask slice number corresponding to T1 images
% function : ROIsMaskWM = fAutoMaskWM(T1map,Tesla,mask)
%
% Authoer : Wanyong Shin
%
% description: WM mask is acquired based on T1 distribution 
% using half-hight cutting off.
%             
% status: stable
%
% versions
%   02-09-2005 (WYS): Initialization
%   05-18-2005 (WYS): updated
%   03-07-2007 (MKS): display match difference

mask = zeros(size(maskSS));
% info1 = dicominfo([path1 '\1.dcm']);
info1 = dicominfo([path1 '\1.dcm'],'Dictionary','dicom-dict.txt'); %Upon the upgrade of the Matlab dicom dictionary
SlicePosition1 = info1.SliceLocation;   
%SlicePosition1 = 18;   %Manually enter LL-EPI slice location
Size1 = double([info1.Rows  info1.Columns]);
PixelSize1 = info1.PixelSpacing';
clear info1

%Grady (12/9/2014)
%N_meas = 200; 

SlicePosition2 = zeros(1,Nslices);
for i = 1:Nslices
%     info2 = dicominfo([path2 '\' num2str(i) '.dcm']);
    info2 = dicominfo([path2 '\' num2str(1+(i-1)*N_meas) '.dcm'],'Dictionary','dicom-dict.txt'); %Upon the upgrade of the Matlab dicom dictionary
    SlicePosition2(i) = info2.SliceLocation;
end

% if (info2.Rows<=112) || (info2.Columns<=112)
%     Size2 = double([info2.Rows  info2.Columns]).*128/112;
%     PixelSize2 = (112/128).*info2.PixelSpacing';
% else
    Size2 = double([info2.Rows  info2.Columns]);
    PixelSize2 = info2.PixelSpacing';
%end

clear info2

% find the corresponding slice position
[a b] = min(abs(SlicePosition2-SlicePosition1));
if a < 0.1
    index = 'match';
    match = 'exact';
%    if ~index
else
    index = 'unmatch';    % changed to call match but unmatch will be noted
    match = 'closest';
    fprintf('Slice location is not matched \n');
    fprintf('Match difference is: %f (T1: %f, DSC: %f)\n',a,SlicePosition1,SlicePosition2(b));
    fprintf('Closest slice number is: %d\n',b);
end
nslice = b;
%nslice = 16;  %gradytest 3/9/2015, slice of LL

% % resize and make WM_DSC mask
if strmatch (index,'match')
	if Size1 == Size2   % if size of matrix is same such as 128x128 <-> 128x128
        if PixelSize1 == PixelSize2
            mask = maskSS;
        else
            mask = imresize(maskSS,2*round([Size1.* PixelSize1./PixelSize2]./2));
            if size(mask,1) > Size1(1)
                cutoffsX = (size(mask,1) - Size1(1))/2;
                mask(size(mask,1)-cutoffsX+1:size(mask,1),:)=[];
                mask(1:cutoffsX,:)=[];
            end
            if size(mask,2) > Size1(2)
                cutoffsY = (size(mask,2) - Size1(2))/2;
                mask(:,size(mask,2)-cutoffsY+1:size(mask,2))=[];
                mask(:,1:cutoffsY)=[];
            end
            fprintf('Different FOV, though, resized...\n');
        end
	else % if 128 x 128 <-> 128 x 104 or 256 x 256
        if PixelSize1 == PixelSize2 
            if (Size1(1) - Size2(1)) > 0
                cutoffsX = (Size1(1) - Size2(1))/2;
                mask = maskSS(cutoffsX+1:size(maskSS,1)-cutoffsX,:);
            elseif (Size1(1) - Size2(1)) < 0
                addX = (Size2(1) - Size1(1))/2;
                mask = [zeros(size(maskSS,1),addY); maskSS; zeros(size(maskSS,1),addY)];
            end
            if (Size1(2) - Size2(2)) > 0
                cutoffY = (Size1(2) - Size2(2))/2;
                mask = maskSS(:,cutoffsY+1:size(maskSS,2)-cutoffsY);
            elseif (Size1(2) - Size2(2)) < 0
                addY = (Size2(2) - Size1(2))/2;
                mask = [zeros(size(maskSS,1),addY) maskSS zeros(size(maskSS,1),addY)];
            end
        elseif Size1.*PixelSize1 == Size2.*PixelSize2
            mask = imresize(maskSS,[Size2]);
        else
            mask = imresize(maskSS,2*round([Size1.* PixelSize1./PixelSize2]./2));
            
            fprintf('error : Different Pixel Size\n');
        end
    end
else
    %fprintf('error : Slice location is not matched')
end

% resize and make WM_DSC mask
if strmatch (index,'match')
    % if FOV is different -> add or delete boundary
    rFOV = (Size1.*PixelSize1)./(Size2.*PixelSize2);
    mask = imresize(maskSS,2*round([Size2.*rFOV]./2));

    if size(mask,1) < Size2(1) 
        addoffsX = (Size2(1) - size(mask,1))/2;
        mask = [zeros(addoffsX,size(mask,2)); mask; zeros(addoffsX,size(mask,2))];
    elseif size(mask,1) > Size2(1) 
        cutoffsX = (size(mask,1) - Size2(1))/2;
        mask = mask(cutoffsX+1:size(mask,1)-cutoffsX,:);
    end
    if size(mask,2) < Size2(2) 
        addoffsY = (Size2(2) - size(mask,2))/2;
        mask = [zeros(size(mask,1),addoffsY) mask zeros(size(mask,1),addoffsY)];
    elseif size(mask,2) > Size2(2) 
        cutoffsY = (size(mask,2) - Size2(2))/2;
        mask = mask(:,cutoffsY+1:size(mask,2)-cutoffsY);
    end
else   
    % if FOV is different -> add or delete boundary
    rFOV = (Size1.*PixelSize1)./(Size2.*PixelSize2);
    mask = imresize(maskSS,2*round([Size2.*rFOV]./2));

    if size(mask,1) < Size2(1) 
        addoffsX = (Size2(1) - size(mask,1))/2;
        mask = [zeros(addoffsX,size(mask,2)); mask; zeros(addoffsX,size(mask,2))];
    elseif size(mask,1) > Size2(1) 
        cutoffsX = (size(mask,1) - Size2(1))/2;
        mask = mask(cutoffsX+1:size(mask,1)-cutoffsX,:);
    end
    if size(mask,2) < Size2(2) 
        addoffsY = (Size2(2) - size(mask,2))/2;
        mask = [zeros(size(mask,1),addoffsY) mask zeros(size(mask,1),addoffsY)];
    elseif size(mask,2) > Size2(2) 
        cutoffsY = (size(mask,2) - Size2(2))/2;
        mask = mask(:,cutoffsY+1:size(mask,2)-cutoffsY);
    end
end

% clear variables
clear path1 path2 maskSS Nslices info1 info2 SlicePosition1 SlicePosition2 Size1  Size2 PixeSize1 PixeSize2 a b