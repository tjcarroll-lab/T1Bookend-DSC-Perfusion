function S_images = fSpatialFilter(images,SF,mask);

% S_images = fSpatialFilter(images,SF,mask);
%
% spatial filtering using convoluntion (smoothing)
% author:   Wanyong Shin
%
% description: image is filtered with SF, which is [Wf, N] 
% smoothing window is N by N
% SF = [Wf,N]
% example N=3 is below 
%                  1           [  1     Wf    1  ]     
%              ----------  x   [  Wf  (Wf^2)  Wf ] 
%             4+4*Wf+Wf^2      [  1     Wf    1  ]  
%
% references
%
% status:   stable 

% versions
%   [04-02-11] (WYS): initial version
%   [05-08-19] (WYS): modified  
%   [07-09-06] (JJM): modified to support any N value in case of Wf = 1

% parse arguments
if nargin < 3
    mask = AutoMask_NS(images);
    if nargin < 2
        SF = [1 1]; 
    end
end

% definition
S_images = images;
Wf=SF(1);   N=SF(2);

if N == 3
    Weigthedwindow = [1 Wf 1;Wf  (Wf^2)  Wf;1 Wf 1]; 
elseif N == 5
    Weigthedwindow = [1 Wf (Wf^2) Wf 1;Wf Wf^2 Wf^3 Wf^2 Wf;Wf^2 Wf^3 Wf^4 Wf^3 Wf^2;Wf Wf^2 Wf^3 Wf^2 Wf;1 Wf (Wf^2) Wf 1];
end

%Jessy (9/6/2007):
if (Wf == 1)
    Weigthedwindow = ones(N,N);
end
%Jessy (end)

%YIJ 20170616 parallel computing
ksize = size(images,3);
istart = round(SF(2)/2);
iend = size(images,1)-round(SF(2)/2)+1;
jstart = round(SF(2)/2);
jend = size(images,2)-round(SF(2)/2)+1;
% if ~matlabpool('size')
%     matlabpool
% end
if SF(2) > 1
    %parfor k = 1:ksize
    for k = 1:ksize
        image = images(:,:,k);      
        for i = istart:iend
            for j = jstart:jend
                if mask(i,j)
                    imagewindow = image(i-(N-1)/2:i+(N-1)/2,j-(N-1)/2:j+(N-1)/2);
                    S_images(i,j,k) = sum(imagewindow(find(isfinite(imagewindow))).*Weigthedwindow(find(isfinite(imagewindow))))/sum(Weigthedwindow(find(isfinite(imagewindow))));
                end
            end
        end
    end
end

    