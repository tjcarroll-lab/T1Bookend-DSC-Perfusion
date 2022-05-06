function [A,alpha,beta,RESNORM,RESIDUAL,EXITFLAG] = gamma_variate_IDC_LeastSquaresFit(times,amplitudes,A,alpha,beta)

%%==================================================================%
%% fitting transit time ditribution curve 
%%==================================================================%
% input:  1) times = x data
%         2) amplitudes = y data
% output: 1) alpha = first gamma variate parameter
%         2) beta = second gamma variate parameter
%         3) RESNORM 
%         4) EXITFLAG
%
% Author: Jessy Mouannes 
% Date:   11-06-2007
% Modified 05-20-2008: changed function to gamma variate IDC model fit
%
%====================================================================%

%% make sure data is in rows
[r,c]=size(times);
if r > c
    times=times';
end; 
[r,c]=size(amplitudes);
if r > c
    amplitudes=amplitudes';
end;

%  if length(amplitudes) == 3
%     newamplitudes = zeros(1,5);
%     newtimes = newamplitudes;
%     newamplitudes(1:2:5)= amplitudes;
%     newamplitudes(2)=mean([amplitudes(1) amplitudes(2)]);
%     newamplitudes(4) = mean([amplitudes(2) amplitudes(3)]);
%     
%     newtimes(1:2:5)= times;
%     newtimes(2)=mean([times(1) times(2)]);
%     newtimes(4) = mean([times(2) times(3)]);
%      
%     amplitudes = newamplitudes; times=newtimes;
%  end
sigmax=find(amplitudes == max(amplitudes)); 
if length(sigmax)>1
    sigmax= sigmax(1);
end
%   if length(amplitudes(1:sigmax)) == 2
%     slope = (amplitudes(sigmax)-amplitudes(1))/(times(sigmax)-times(1));
%     addsignal = zeros(1,3);
%     addtime = addsignal;
%     addsignal(1:2:3)= amplitudes(1:sigmax);
%     addsignal(2)=0.8*(amplitudes(sigmax)-amplitudes(1))+amplitudes(1);
%     
%     
%     addtime(1:2:3)= times(1:sigmax);
%     addtime(2)=(addsignal(2)-amplitudes(1))/slope;
%     
%      
%     amplitudes = [addsignal amplitudes(sigmax+1:end)]; 
%     times = [addtime times(sigmax+1:end)];
%  end
 


%%% fit function %%%
fun = inline('x(3).*xdata.^x(1).*exp(-xdata./x(2))','x','xdata');
% fun = inline('1./(x(2).^x(1).*gamma(x(1)).*xdata.^(x(1)-1).*exp(-xdata./x(2)))','x','xdata');

%%%estimate
% X0 = [alpha,beta,A,shift];
% % lowB = [0,0];
% % highB = [5*1.5,5*1.8];
% lowB = [0,0,0,0];
% highB = [10,10,10,2];
% %%%perform fit
% pin = [A alpha beta];
X0 = [alpha,beta,A];
lowB = [0,0,0];
highB = [10,10,10];

if strcmp(version('-release'),'2013a') 
    options = optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt','MaxFunEvals',15000);
elseif  strcmp(version('-release'),'2012b')
    options = optimset('Display','on','Algorithm','levenberg-marquardt','MaxFunEvals',15000);
else
    %options = optimset('Display','on','LevenbergMarquardt','on');
    %options = optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt','MaxFunEvals',15000);
    options = optimset('Display','off','Algorithm','trust-region-reflective','MaxFunEvals',15000);
end
% err = ones(length(amplitudes),1);
% fcp(1) = 0.01; fcp(2) = 20; fcp(3) = 0.001; 
% notfixed = [1;1;1];
% 
% 
% [p,sig,f,covp,corp,r2,rv,interp_time] = levmarq_lsqr_fit_PV(times,amplitudes,err,pin,notfixed,'gamma_var_fun',fcp);
% 
% %try
% alpha = p(2);
% beta = p(3);
% strength = p(1);
% shift = 0;
% f = feval('gamma_var_fun',times,p);
% 
% error.percent_error = abs((f(1:end)-amplitudes(1:end))./amplitudes(1:end));
% %catch
% %    1;
% %end
% error.percent_error = sum((error.percent_error(1:end)))/length(error.percent_error(1:end));
% error.sum_error = abs((f(1:end)-amplitudes(1:end)));
% error.sum_error = mean(error.sum_error);

[X,RESNORM,RESIDUAL,EXITFLAG]=lsqcurvefit(fun,X0,times,amplitudes,lowB,highB,options);
A = X(3);
alpha = X(1);
beta = X(2);

%%% write values %%%
%if EXITFLAG
    
%else
%    alpha = 0;
%    beta = 0;
%    strength = 0;
%    shift = 0;
%end

% clear variabels
clear times amplitudes fun X0 X 