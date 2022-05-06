function CF = fWaterExCorr_MS2(R1,Tesla)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Water exchange correction based on MS cases (Toronto)
%
% function CF = fWaterExCorr_MS(R1,Tesla)
%
% Model determined by fitting experimental data to 
% to a second degree polynomial function, similar to 
% the Hazlewood Model - This is using the CBV vs. dR1 
% data from the NAWM ROI selected for calibration 
%
% Author: Jessy Mouannes
% Date: August 31, 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global injectionNum;

if  strcmp(Tesla,'3.0T')
    %CF = 2.2334e-014*R1^2 + 0.2539*R1 + 0.0792;
    %CF =  4.1815e-014*R1.^2 + 0.3419*R1 + 0.0792; %grady test
    %CF = .01675*R1.^2 + .4035*R1 + .015;
    switch injectionNum
        case 1
            disp(['WCF Injection number1: ' num2str(injectionNum)]);
            %CF = 0.3421*R1.^2 + 0.1572*R1 + 7.827e-12; %yong
            %CF = 0.09175*R1.^2 + 0.0335*R1 + 0.1; %grady volunteer
            %CF = 0.4415*R1.^2 + 9.932e-14*R1 + 9.548e-14;
            %CF = 0.3757*R1.^2 + 0.1363*R1 + 5.313e-09; %20171005
            %CF = 0.2051*R1.^2 + 0.3569*R1 + 1.942e-08; %using 2.2 from water ex simulation
            
            CF = 0.3539.*R1.^2 + 0.1106.*R1 + 0.0601; %using 0.02 from sim - 20171031 stable for paper
            %CF = 4.36e-12.*R1.^2 + 0.3691.*R1 + 0.2499; % YIJ 20180718 new fit experimental
            
            %CF = (0.2601.*R1.^2 + 0.1977.*R1 + 0.1187); %using 0.022 from sim
            %CF = (0.403.*R1.^2 + 0.01059.*R1 + 0.08402); %using 0.02 from sim 20171101
            %CF = (0.1428.*R1.^2 + 0.4585.*R1 + 0.0003805); %using 0.022 from sim 20171101 (3,4,5,6,7)
        case 2
            disp(['WCF Injection number2: ' num2str(injectionNum)]);
            %CF = 0.17075*R1.^2 + 0.00001*R1 + 0; %grady volunteer
            %CF = .02775*R1.^2 + .5935*R1 + .01; %grady patient
            %CF = 0.124*R1.^2 + 1.544e-09*R1 + 0.2458; %yong using '2'
            %CF = 0.31*R1.^2 + 1.404e-08*R1 + 0.6144; %yong using '5'
            %CF = 0.007004*R1.^2 + 0.1468*R1 + 0.05783;
            %CF = 0.03492*R1.^2 + 3.62e-09*R1 + 0.2517; %20171005
            %CF = 0.1013*R1.^2 + 2.178e-11*R1 + 0.73; %using 5.8 from water ex simulation
            CF = (0.2174.*R1.^2 + 1.994e-06.*R1 + 0.4875);%.*0.55; %using 0.036 from sim - 20171030 stable for paper (assuming carbogen)
            %CF = (0.02511.*R1.^2 + 0.5282.*R1 + 0.2475); %using 0.04 from sim
            %CF = (2.221e-14.*R1.^2 + 0.3264.*R1 + 0.05355); %using 2
            %CF = (0.00736.*R1.^2 + 0.7348.*R1 + 0.02409); %using 0.042 from sim 20171101
            %CF = (0.1365.*R1.^2 + 0.2208.*R1 + 0.5502); %using 0.044 from sim 20171101 (3,4,5,6,7)
        otherwise
            disp(['WCF Injection number ("otherwise"): ' num2str(injectionNum)]);
            %CF = 0;
            CF = 0.3539.*R1.^2 + 0.1106.*R1 + 0.0601;
            %warning(['No WCF determined for injection number ' num2str(injectionNum)]);
    end               
else
    CF = 0;
end