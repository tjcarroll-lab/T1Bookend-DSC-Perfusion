function T = TimeTableforIRLLEPI(path,Length_vec)

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
header = dicominfo([path '\1.dcm'],'Dictionary','dicom-dict.txt'); %Upon the upgrade of the dicom dictionary


% useful time parameter to calculate normal and corrected time table
TG = 7.747+header.EchoTime;
%TG = header.EchoTime;
TR = header.RepetitionTime;
%TR = 6*header.EchoTime;
%TR = (73/120)*1000;
%TR = 45.3697;
%TR = 45.3697+header.EchoTime;

T  = TG:TR:TG+TR*(Length_vec-1);

clear TR TE header