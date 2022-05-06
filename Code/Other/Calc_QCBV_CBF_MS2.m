function [CF,qimages] = Calc_QCBV_CBF_MS2(images,T1s,fitting,Tesla);


%                                                               02-11-2004 
global injectionNum;
global glblTargetPath;
global seqType;

% physiological factor
CF.PhysioK = 55/75/1.04; % H_LV = 0.45 H_SV = 0.25

% make temp mask to avoid NaN and Inf
R1 = 1000*(1/T1s.blood_post-1/T1s.blood_pre);

%Jessy (3/31/2009): MS cases water exchange correction model
CF.WCF = fWaterExCorr_MS2(R1,Tesla);
% CF.WCF = fWaterExCorr(R1,1/1.1,'WM','fast',Tesla);


CF.ratio_CBV = fitting.CBVssFast.b1/fitting.CBVdsc.b1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%grady change
% if (fitting.CBVssFast.b1 > 1.8)
%     CF.ratio_CBV = 1.8/fitting.CBVdsc.b1;
% end
% if (fitting.CBVssFast.b1 < .5)
%     CF.ratio_CBV = .5/fitting.CBVdsc.b1;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch injectionNum
    case 1
        disp(['CBV Injection number1: ' num2str(injectionNum)]);
        CF.CBVsswm_ratio = 1;
    case 2
        disp(['CBV Injection number2: ' num2str(injectionNum)]);
        tmp = load([glblTargetPath '\Result_MSwcf2\P001' seqType '.mat']);
        %CF.CBVsswm_ratio = fitting.CBVssFast.b1 / tmp.ROIs.data.fitting.CBVssFast.b1;
        %CF.CBVsswm_ratio = fitting.Conc.b1 / tmp.ROIs.data.fitting.Conc.b1;
        %CF.CBVsswm_ratio = (fitting.CBVssFast.b1 / tmp.ROIs.data.fitting.CBVssFast.b1) * (fitting.Conc.b1 / tmp.ROIs.data.fitting.Conc.b1);
        CF.CBVsswm_ratio = 1;
    otherwise
        disp(['CBV Injection number ("otherwise"): ' num2str(injectionNum)]);
        %CF.CBVsswm_ratio = 0;
        CF.CBVsswm_ratio = 1;
        warning(['Injectio number is greater than 2. injection number ' num2str(injectionNum)]);
end
    
qimages = CF.PhysioK*CF.WCF*CF.ratio_CBV*images*CF.CBVsswm_ratio;


% clear variables
clear images CBVss CBVdsc Tesla R1 