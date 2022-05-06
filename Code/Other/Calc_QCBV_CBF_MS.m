function [CF,qimages] = Calc_QCBV_CBF_MS(images,T1s,fitting,Tesla);


%                                                               02-11-2004 

% physiological factor
CF.PhysioK = 55/75/1.04;

% make temp mask to avoid NaN and Inf
R1 = 1000*(1/T1s.blood_post-1/T1s.blood_pre);

%Jessy (3/31/2009): MS cases water exchange correction model
CF.WCF = fWaterExCorr_MS(R1,Tesla);
% CF.WCF = fWaterExCorr(R1,1/1.1,'WM','fast',Tesla);



CF.ratio_CBV = fitting.CBVssFast.b1/fitting.CBVdsc.b1;

qimages = CF.PhysioK*CF.WCF*CF.ratio_CBV*images;


% clear variables
clear images CBVss CBVdsc Tesla R1 