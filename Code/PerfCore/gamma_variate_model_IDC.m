function [contrast] = gamma_variate_model_IDC(t,strength,alpha,beta)
%-----------------------------------------------------------%
% gamma_variate:                                            %
% plot contrast = (strength)*(t**alpha)*exp(-beta*t)        %
%                                                           %
%                          Timothy J. Carroll               %
%                          March 15, 2002                   %
%                                                           %
% Modified by: Jessy Mouannes                               %
% Date: May 20, 2008                                        %
% Modification: gamma variate model for indicator dilution  %
% curve (IDC) fitting                                       %
%-----------------------------------------------------------%

% contrast = strength.*(t.^alpha).*exp(-1.0*beta*t);
shift = 0;
contrast = strength.*((t+shift).^alpha).*exp(-1.0./beta.*(t+shift));
%
% the next few lines are usefull for debugging;
%contrast = t.^alpha;
%contrast = exp(-1.0*beta*t);
%function [contrast] = gamma_variate(t,strength,alpha,beta)
