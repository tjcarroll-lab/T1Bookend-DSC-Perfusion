function [T1,M0,InvF,RESNORM,signalFit] = fIR_LLEPI_fitting(times,amplitudes)

%%========================================================================%
%%      fitting tfiti curve up to | A - B*exp(-t/tau) |  
%%========================================================================%
% input :   1) T (vector)         : time tables    
%           2) Intensity (vector) : tfiti data curve
%%========================================================================%
% output:   1) tau (scalar)     : T1
%           2) error (scalar)   : error in tau
%           3) M0 (scalar)      : initial magnization magnitude
%                                                                02-17-2004
                                                             

% make sure data is in rows
[r,c]=size(times);
if r > c
    times=times';
end; 
[r,c]=size(amplitudes);
if r > c
    amplitudes=amplitudes';
end;

%%% fit function %%%
fun=inline('abs(x(1) - x(2)*exp(-xdata/x(3)))','x','xdata');

%%% estimate M0 %%%
M00=mean(amplitudes(1)); 

%%% estimate T1 %%%
[A,IX]=min(amplitudes);  % get null point amplitude and index
T10=1.44*times(IX); % estimate T1

%%% initial fit parameters %%%
X0=[M00,2*M00,T10]; 
%lowB = [M00-1,10,20];
%highB = [M00+1,12000,3000];
lowB = []; highB = [];
%%% perform fit %%%
options = optimset('Display','off','Algorithm','levenberg-marquardt','MaxFunEvals',15000);
[X,RESNORM,RESIDUAL,EXITFLAG]=lsqcurvefit(fun,X0,times,amplitudes,lowB,highB,options);
%EXITFLAG = 0; RESNORM = 0; RESIDUAL = 0;
%%% write values %%%
if EXITFLAG
    tau=X(3);    
    M0=X(1);
    %M0=X(2)-X(1);
    InvF = X(2)/X(1);
    T1 = tau*(InvF-1);
    %solution for improper sequence from 2009 for AvivFinalNov2014
    %T1 = real(-X(3)*log(.37*X(1)/X(2)));
    if  T1 < 0
        T1 = 0;
    end
    if T1 > 2000
        T1 = 2000;
    end
    for i = 1:size(times,2)
        signalFit(i) = abs(X(1) - X(2)*exp(-times(i)/X(3)));
    end
else
    M0=0;    InvF = 0;    T1 = 0; signalFit = amplitudes;
end
% clear variabels
clear times amplitudes fun M00 X0 A IX X RESIDUAL