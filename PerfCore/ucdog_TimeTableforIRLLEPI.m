function T = ucdog_TimeTableforIRLLEPI(path,Length_vec)
% adapted from TimeTableforIRLLEPI
% times are stored as TriggerTime in philips scanner
%
% Author: Yong Ik Jeong
% Date: 2016-03-22
%%%%===================================================================%%%%
%% make a time tables for IR LL-EPI sequence
%%%%===================================================================%%%%
%%  Input   1) path    : Intensity vector
%%          2) DT   (structure) : info by reading dcm header file
%%  Output  1) T        (vector)    : Time table 
%%%%===================================================================%%%%
%%              Created by Wanyong Shin                                  %%
%%                                  March 16 2005                        %%
%%%%===================================================================%%%%

% read info
% header = dicominfo([path '\1.dcm']);
T = [];
for ii = 1:Length_vec
    header = dicominfo(fullfile(path,[num2str(ii) '.dcm']),'Dictionary','dicom-dict.txt'); %Upon the upgrade of the dicom dictionary
    T = [T header.TriggerTime];
end

% % YIJ 20180807 test
% T = linspace(26,2066,256);
% T
% useful time parameter to calculate normal and corrected time table
% TG = header.TriggerTime;
% TR = header.RepetitionTime;
% T  = TG:TR:TG+TR*(Length_vec-1);

clear TR TE header