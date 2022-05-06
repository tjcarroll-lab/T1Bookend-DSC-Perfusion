function AIF_struct = Calc_AIFstruct(path_DSC,AIF,n_slice_AIF,N_slices)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function AIF_struct = Calc_AIFstruct(AIF,n_slice_AIF,N_slices)
% 
% This function calculate a shifted AIF version coresponding to each
% slice's order of acquisition
%
% Input:  1) AIF: struct containing the fitted AIF parameters
%         2) n_slice_AIF: number of slice containing the AIF voxels
%         3) N_slices: number of acquired slices
%
% Output: 1) AIF_struct: struct of the size of the number of slices
%
% Author: Jessy Mouannes
% Date: June 26, 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dcmfiles = dir([path_DSC '/*.dcm']);
% prevslc = 0;
% for ii = 1:length(dcmfiles)
%     tmphdr = dicominfo([path_DSC '/' num2str(ii) '.dcm']);
%     hdrs{ii} = tmphdr;
%     %currslc = tmphdr
%     tmphdr.InstanceNumber
%     tmphdr.SliceLocation
%     tmphdr.TemporalPositionIdentifier
%     tmphdr.AcquisitionTime
%     tmphdr.InstanceCreationTime
%     fprintf('------------------------------\n');
% end

%% Compute extended AIF concentration vector by length = TR
header = dicominfo([path_DSC '/1.dcm']); 
Dt = (header.RepetitionTime)/1000; % from ms to sec
time = 0:Dt/N_slices:(AIF.RTP-AIF.BATP+1)*Dt; % +1 to account for the additional TR due to shift 
AIF_CtExt = gamma_variate_model_IDC(time,AIF.FitParam.A,AIF.FitParam.alpha,AIF.FitParam.beta);

%% Compute AIF_struct

TE = header.EchoTime;

%Index vector for re-adjusting AIF size
AIF_length = AIF.RTP-AIF.BATP+1;
Index_vec = ones(1,AIF_length);
for i = 2:AIF_length
    Index_vec(i) = N_slices*(i-1)+1;
end

%% YIJ 20170913 Assume ascending ordering 
AIF_temp = AIF.Ct;
AIF_struct{n_slice_AIF} = AIF_temp(Index_vec);%*TE;
n_slice = n_slice_AIF + 1;
while (n_slice <= N_slices)
    AIFshift = n_slice - n_slice_AIF;
    AIF_temp = AIF_CtExt(1+AIFshift:(length(AIF.Ct)+AIFshift));
    AIF_struct{n_slice} = AIF_temp(Index_vec);%*TE;
    n_slice = n_slice + 1;
end

n_slice = n_slice_AIF - 1;
while n_slice
    AIFshift = n_slice - n_slice_AIF;
    AIF_temp = [zeros(1,abs(AIFshift)) AIF_CtExt(1:(length(AIF.Ct)+AIFshift))];
    AIF_struct{n_slice} = AIF_temp(Index_vec);%*TE;
    n_slice = n_slice - 1;
end

% %% Assume some interleaved ordering
% % Even slice number
% if (mod(n_slice_AIF,2)==0) 
%     %AIF slice
%     AIF_temp = AIF.Ct;
%     AIF_struct{n_slice_AIF} = AIF_temp(Index_vec)*TE;
%     %Incrementing over even slices
%     n_slice = n_slice_AIF + 2;
%     while (n_slice <= N_slices)
%         AIFshift = (n_slice-n_slice_AIF)/2;
%         AIF_temp = AIF_CtExt(1+AIFshift:(length(AIF.Ct)+AIFshift));
%         AIF_struct{n_slice} = AIF_temp(Index_vec)*TE;
%         n_slice = n_slice + 2;
%     end
%     %Decrementing over even slices
%     n_slice = n_slice_AIF - 2;
%     while n_slice
%         AIFshift = (n_slice-n_slice_AIF)/2;
%         AIF_temp = [zeros(1,abs(AIFshift)) AIF_CtExt(1:(length(AIF.Ct)+AIFshift))];
%         AIF_struct{n_slice} = AIF_temp(Index_vec)*TE;
%         n_slice = n_slice - 2;
%     end
%     %Decrementing over all odd slices
%     if (mod(N_slices,2)==0) %even number of slices
%         n_slice = N_slices - 1;
%     else %odd number of slices
%         n_slice = N_slices;
%     end
%     AIFskip = (n_slice_AIF-2)/2;
%     AIFshift = 1;
%     while (n_slice >= 1)    
%         AIF_temp = [zeros(1,(AIFshift+AIFskip)) AIF_CtExt(1:(length(AIF.Ct)-AIFshift-AIFskip))];
%         AIF_struct{n_slice} = AIF_temp(Index_vec)*TE;
%         n_slice = n_slice - 2;
%         AIFshift = AIFshift + 1;
%     end
%     
% %Odd slice number   
% else 
%     %AIF slice
%     AIF_temp = AIF.Ct;
%     AIF_struct{n_slice_AIF} = AIF_temp(Index_vec)*TE;
%     %Decrementing over odd slices
%     n_slice = n_slice_AIF - 2;
%     while (n_slice >= 1)
%         AIFshift = (n_slice-n_slice_AIF)/2;
%         AIF_temp = [zeros(1,abs(AIFshift)) AIF_CtExt(1:(length(AIF.Ct)+AIFshift))];
%         AIF_struct{n_slice} = AIF_temp(Index_vec)*TE;
%         n_slice = n_slice - 2;
%     end
%     %Incrementing over odd slices
%     n_slice = n_slice_AIF + 2;
%     while (n_slice <= N_slices)
%         AIFshift = (n_slice-n_slice_AIF)/2;
%         AIF_temp = AIF_CtExt(1+AIFshift:(length(AIF.Ct)+AIFshift));
%         AIF_struct{n_slice} = AIF_temp(Index_vec)*TE;
%         n_slice = n_slice + 2;
%     end
%     %Incrementing over all even slices
%     if (mod(N_slices,2)==0) %even number of slices
%         AIFskip = (N_slices-n_slice_AIF-1)/2;
%     else %odd number of slices
%         AIFskip = (N_slices-n_slice_AIF)/2;
%     end
%     n_slice = 2;
%     AIFshift = 1;
%     while (n_slice <= N_slices)
%         AIF_temp = [AIF_CtExt((1+AIFshift+AIFskip):(length(AIF.Ct)+AIFshift+AIFskip))];
%         AIF_struct{n_slice} = AIF_temp(Index_vec)*TE;
%         n_slice = n_slice + 2;
%         AIFshift = AIFshift + 1;
%     end
% end

