function [] = debug_perf(dsc, header, voxelloc, result)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dsc = double(dsc);
Dt = (header.RepetitionTime)/1000; % from ms to sec
echoTime = header.EchoTime;
Wij_SVD = 0.2;
Wij = 0.2;
[path_DSC,~,~] = fileparts(header.Filename);
AIFslice = result.ROIs.data.AIFslice;
n_slice_AIF = result.ROIs.positions.n_slice_AIF;
[nr, nc, N_meas, N_slices] = size(dsc);

sampdelayfix = 1;
DDfix = 1;

AIF_struct = Calc_AIFstruct(path_DSC,AIFslice,n_slice_AIF,N_slices);
AIF = AIFslice;

CmtWM = result.cutoffs_ROIs.CmtWM;

VEIN = result.ROIs.data.VEIN;
ATDmap = result.images{strcmpi(result.image_names, 'ATDmap')};
IBATmap = result.images{strcmpi(result.image_names, 'IBATmap')};
BATmap = result.images{strcmpi(result.image_names, 'BATmap')};
BRTmap = result.images{strcmpi(result.image_names, 'BRTmap')};

k = voxelloc(3);

DDfit = [];

AIF.Ct = AIF_struct{k};

% get delay and dispersion fit using resampled AIF
[DDfit] = DDcorr_test(AIF,Dt,Wij_SVD,VEIN);


images = squeeze(dsc(:,:,:,k));

%YIJ 20170615 smooth mask
mask = automaskns(imfilter(images,fspecial('gaussian',[3 3],2)));

% spatial filtering
S_images = fSpatialFilter(images,[2 3],mask);
% 20190725 YIJ: Try filter twice (might be better for Tmax)
%S_images = fSpatialFilter(S_images,[2 3],mask);


% define variables
Kh = (1-0.45)/(1-0.25);
psi = 1.04; % mg/ml  --> g/mL?

if isfield(AIF,'ITP')
    ITP=AIF.ITP; BAT=AIF.BATP; RTP=AIF.RTP;
else
    ITP=2; BAT=AIF.BATP; RTP=AIF.RTP;
end

warning off;

% save original AIF
oAIF = AIF;

i = voxelloc(1);
j = voxelloc(2);

% call info. from Cmt strucure
Signal = squeeze(S_images(i,j,:))';

isbolus = boluscheck(Signal,oAIF,[IBATmap(i,j,k) BATmap(i,j,k) BRTmap(i,j,k)], echoTime);

%YIJ 20170613
Signal = smooth(Signal,5,'sgolay',3)';

if DDfix && BATmap(i,j,k) && IBATmap(i,j,k) && BRTmap(i,j,k)
    % get regional conc. curve
    Cmt = Sig2Conct_AIF(Signal, [IBATmap(i,j,k) BATmap(i,j,k) BRTmap(i,j,k)], echoTime);
    Cmt2 = Sig2Conct_AIF(Signal, [IBATmap(i,j,k) BAT BRTmap(i,j,k)], echoTime);
    % apply DDfix to AIF
    AIF = applyDDfix_test(oAIF,Dt,ATDmap(i,j,k),DDfit,Wij);
    %tmp = AIF;
    AIF.Ct = AIF.Ct .* sum(oAIF.Ct)/sum(AIF.Ct);
    if abs(sum(oAIF.Ct) - sum(AIF.Ct)) > .0000001
        1;
    end
    if ATDmap(i,j,k) > 5
        1;
    end
else
    %Cmt = Sig2Conct_AIF(Signal, [CmtWM(1:3)], echoTime);
    Cmt = Sig2Conct_AIF(Signal, [IBATmap(i,j,k) BAT BRTmap(i,j,k)], echoTime);
    Cmt2 = Sig2Conct_AIF(Signal, [IBATmap(i,j,k) BAT 200], echoTime);
    %                 Signal = smooth(Signal,.05,'sgolay',1)';
    %                 Cmt2 = Sig2Conct_AIF(Signal, [IBATmap(i,j) BAT BRTmap(i,j)], echoTime);
    AIF = oAIF;
end

% calculate Residue function
fResidue = fResidue_test(AIF,Cmt,Dt,Wij);
fResidue2 = fResidue_test(AIF,Cmt2,Dt,Wij);

[MaxIn MaxN] = max(fResidue);
[MaxIn2 MaxN2] = max(fResidue2);

% clear variables
clear images mask
% clear variables
clear header Dt Wij_SVD AIF mask images Wf S_images temp_CBV_DSC temp_CBF_SVD temp_CMTT_CVP tempTmaxMap temp_dBATmap

end

function Residue = fResidue_test(AIF,Cmt,Dt,Wij)

% in case that AIF is not structure
if isfield(AIF,'Ct')
    temp_AIF.Ct = AIF.Ct;
else
    temp_AIF.Ct = AIF;
end
if isfield(AIF,'BATP')
    temp_AIF.BATP = AIF.BATP;
else
    temp_AIF.BATP = 1;
end
if isfield(AIF,'RTP')
    temp_AIF.RTP = AIF.RTP;
else
    temp_AIF.RTP = length(temp_AIF.Ct);
end

% smoothing AIF 
N = length(temp_AIF.Ct);

%input_AIF = input_AIF(temp_AIF.BATP:temp_AIF.RTP);
windowdiff = length(Cmt.BATP:Cmt.RTP) - length(temp_AIF.BATP:temp_AIF.RTP);
if windowdiff < 0
    %error('window diff negative');
    addaif = 0;
    addct = abs(windowdiff);
else
    addaif = windowdiff;
    addct = 0;
end

if windowdiff ~= 0
    1;
end

input_AIF = temp_AIF.Ct;
input_AIF = [input_AIF zeros(1,addaif)];

if length(Cmt.BATP:Cmt.RTP) == length(Cmt.Ct)
    input_Cmt = [Cmt.Ct zeros(1,addct)];
else
    input_Cmt = [Cmt.Ct(Cmt.BATP:Cmt.RTP) zeros(1,addct)];
end

n = length(input_AIF);

% make AIF matrix
for i = 1:n
    for j=1:i
        AIF_Mat(i,j) = Dt*input_AIF(i-j+1) ;
    end
end 

% sigular value decomposition
[U,S,V] = svd(AIF_Mat);
inv_S = invert_S(S,Wij); 
invA_mat = V*inv_S*U';

% calculation Residue function
Residue = (invA_mat*input_Cmt')';

figure;plot(input_AIF);hold on;plot(input_Cmt);
figure;plot(Residue);

% clear variables
clear AIF Cmt Dt Wij input_AIF input_Cmt AIF_Mat U S V inv_S invA_mat 

end
