%% Example code for sorting DSC and LL dicoms into correct format and order
addpath('.\Other');

sortPerf_ICAD_UC('D:\Users\CarrollLab\Desktop\Yongs Research\Example Cases\Perfusion\P001\ep2d_perf');
sortPerf_ICAD_UC('D:\Users\CarrollLab\Desktop\Yongs Research\Example Cases\Perfusion\P001\LL_EPI_POST');
sortPerf_ICAD_UC('D:\Users\CarrollLab\Desktop\Yongs Research\Example Cases\Perfusion\P001\LL_EPI_PRE');
%% Example code for loading and viewing DSC dicoms (for AIF/VOF ROI)
addpath('.\Other');

dsc = loadDSC('D:\Users\CarrollLab\Desktop\Yongs Research\Example Cases\Perfusion\P001\ep2d_perf');
AIFviewer(dsc);
%% Example code for running perfusion code
addpath('.\PerfCore');

Auto_qCBF_Philips_DD_all('D:\Users\CarrollLab\Desktop\Yongs Research\Example Cases\Perfusion');