# Perfusion instructions 

This code takes repositories and returns both delay and dispersion correct and non delay and dispersion corrected CBF maps.

1. The perfusion code expects the perfusion dicoms to be in a specific directory structure as follows:
   - Case1/ (the name of this folder is arbitrary)
      - P001/ (the code looks for P### folders)
          - ep2d_perf/ (contains the DSC dicoms)
          - LL_EPI_POST/ (contains the T1 Look-Locker post-contrast)
          - LL_EPI_PRE/ (contains the T1 Look-Locker pre-contrast)

2. All the dicoms (DSC and the LLs) need to be renamed into 1.dcm, 2.dcm, 3.dcm, ... format. Also DSC dicoms need to reordered such that the dicoms are sorted by time first then slice (e.g. s0t0, s0t1, s0t2, ..., s1t0, s1t1, s1t2, ...). Use sortPerf_ICAD_UC.m function to rename and sort the dicoms.

3. Before running perfusion code, AIF and VOF should be manually chosen for best results (see 'Manual AIF & VOF instructions.txt'). But the code can run without it (though it's possible the automated AIF/VOF will fail). Run "Auto_qCBF_Philips_DD_all('path to directory')".
   - 'path to directory': this is the path to the directory containing the 'P001' folder

4. When it finishes, perfusion results are stored in the 'Result_MSwcf2' folder (if using T1 quantification), or in the 'DSCanalysis' folder (if relative) as a .mat file.

* See run_perf_script.m for example

# Manual AIF & VOF Instructions


1. Load the DSC dicoms by using "dsc = loadDSC('folder containing DSC')".

2. Use 'AIFviewer.m' function to browse through DSC ("AIFviewer(dsc)").

3. Click on the 'ROI AIF' button to draw an ROI, then click on 'Save ROI'.

4. ROI will be saved in MATLAB workspace as 'aif_roi'. Rename this variable to 'aifmask' and save it as 'AIF_Mask_P001GE_M.mat' in the folder that contains the 'P001' folder.

5. Draw an ROI for VOF and rename the new 'aif_roi' to 'veinmask'. And save it as 'Vein_Mask_P001GE_M.mat' in the same folder as above.
