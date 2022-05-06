function [S_out] = invert_S(S_in,wij)

%%========================================================================%
%%       Invert the diagnonal SVD S matrix and apply the wij theshold 
%%========================================================================%
% input :   1) S_in (ixj matrix)    : singular matrix from SVD
%           2) Wij (scalar)         : Threshold in SVD
%%========================================================================%
% output:   1) S_out (vector)       : 1./S_in in case that S_in > Wij*S_max
%                                                               02-19-2004 

if nargin < 2
    wij = 0;
end

% useful variables
S_max = max(max(S_in));

for i = 1:size(S_in,1)
   if S_in(i,i) > wij.*S_max
      S_out(i,i)= 1./S_in(i,i);
   else
      S_out(i,i) = 0.0;
   end
end

clear S_in wij S_max