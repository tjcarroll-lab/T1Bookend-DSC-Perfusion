function [best_fit] = gamma_variate_IDC_gridsearch(time, signal);
%---------------------------------------------------------------------------------%
% function [best_fit] = gridsearch(time, signal);
%
% Simple brute force grid search fitting algorithm. Note fast or elegant, but     %
% it will help validate some of the other fitting.                                %
%                                                                                 %
%    best_fit(1) = sum of squared differences
%    best_fit(2) = A (Amplitude)                                                  %
%    best_fit(3) = alpha 
%    best_fit(4) = beta 
%
%    Author: Jessy Mouannes
%    Date: May 20, 2008
%    modified by parmede to improve results for fitting AIF to small datasets   %
%---------------------------------------------------------------------------------%
sigmax=find(signal == max(signal)); 
if length(sigmax)>1
    sigmax = sigmax(1);
end
%  if length(signal(1:sigmax)) == 2
%     slope = (signal(sigmax)-signal(1))/(time(sigmax)-time(1));
%     addsignal = zeros(1,3);
%     addtime = addsignal;
%     addsignal(1:2:3)= signal(1:sigmax);
%     addsignal(2)=0.8*(signal(sigmax)-signal(1))+signal(1);
%     
%     
%     addtime(1:2:3)= time(1:sigmax);
%     addtime(2)=(addsignal(2)-signal(1))/slope;
%     
%      
%     signal = [addsignal signal(sigmax+1:end)]; 
%     time = [addtime time(sigmax+1:end)];
%  end
%  


  best_fit    = zeros(4,1);
  best_fit(1) = 99999999.9;
  grid_vector = zeros(1,size(signal,2));

%Iterations  
%alpha_step = 0.0000001;
alpha_vect = 0:0.1:10;
beta_vect = 0:0.1:10;
A_vect = 0:0.005:0.05;
% alpha_vect = 0:1:10;%[0:0.1:10];
% beta_vect = 0:1:10;%0:0.1:10;
% A_vect = 0:0.005:0.05;%0:.005:0.05;
shift_vect = 0;

grid_size = size(alpha_vect,2)*size(beta_vect,2)*size(A_vect,2);
count = 0;

for i = 1:size(alpha_vect,2)
    alpha = alpha_vect(i);
    for j = 1:size(beta_vect,2)
        beta = beta_vect(j);
        for k = 1:size(A_vect,2)
            A = A_vect(k);
            for l = 1:size(shift_vect,2)
                shift = shift_vect(l);
                count = count + 1;
                grid_vector = gamma_variate_model_IDC(time,A,alpha,beta);
                sum_sq      = sum((signal-grid_vector).^2);
                %             As(count) = A;
                %             alphas(count) = alpha;
                %             betas(count) = beta;
                if sum_sq < best_fit(1)
                    best_fit(1) = sum_sq;
                    best_fit(2) = A;
                    best_fit(3) = alpha;
                    best_fit(4) = beta;
                    best_fit(5) = shift;
                end
            end
        end;
    end;
end;
% [a b] = min(sum_sq);
% best_fit(1) = sum_sq(b);
% best_fit(2) = As(b);
% best_fit(3) = alphas(b);
% best_fit(4) = betas(b);
count;
