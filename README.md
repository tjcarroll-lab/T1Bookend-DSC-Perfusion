# DSC Perfusion instructions (mac)

This code takes repositories and returns both delay and dispersion correct and non delay and dispersion corrected CBF maps for a human scan. It has been used in papers including: 

1) Jeong YI, Christoforidis GA, Saadat N, Kawaji K, Cantrell CG, Roth S, Niekrasz M, Carroll TJ. Absolute quantitative MR perfusion and comparison against stable-isotope microspheres. Magn Reson Med. 2019 Jun;81(6):3567-3577. doi: 10.1002/mrm.27669. Epub 2019 Feb 8. PMID: 30737833; PMCID: PMC6435411.

Note: for a new session in matlab, remember to addpath. 
e.g. 
>> addpath '/Users/neuroimaging/Desktop/MR-Code/DSC_code/Perfusion/Code/PerfCore'\
>> addpath '/Users/neuroimaging/Desktop/MR-Code/DSC_code/Perfusion/Code/Other'


## Pre-Processing preparatory requirements
1. The perfusion code expects the perfusion dicoms to be in a specific directory structure as follows:
   - Case1/ (the name of this folder is arbitrary)
      - P001/ (the code looks for P### folders)
          - ep2d_perf/ (contains the DSC dicoms)
          - LL_EPI_POST/ (contains the T1 Look-Locker post-contrast)
          - LL_EPI_PRE/ (contains the T1 Look-Locker pre-contrast)

   
2. All the dicoms (DSC and the LLs) need to be renamed into 1.dcm, 2.dcm, 3.dcm, ... format. Also DSC dicoms need to reordered such that the dicoms are sorted by time first then slice (e.g. s0t0, s0t1, s0t2, ..., s1t0, s1t1, s1t2, ...). 
   Use sortPerf_ICAD_UC.m function to rename and sort the dicoms.
   Or use the python notebook DSCPerfusionSorter_WIP.ipynb. 

3. Before running perfusion code, AIF and VOF should be manually chosen for best results (see 'Manual AIF & VOF instructions.txt'). But the code can run without it (though it's possible the automated AIF/VOF will fail). 

4. When it finishes, perfusion results are stored in the 'Result_MSwcf2' folder (if using T1 quantification), or in the 'DSCanalysis' folder (if relative perfusion) as a .mat file.

* See run_perf_script.m for example

## Manual AIF & VOF Instructions

1. Load the DSC dicoms by using "dsc = loadDSC('folder containing DSC')".

2. Use 'AIFviewer.m' function to browse through DSC ("AIFviewer(dsc)").

3. Click on the 'ROI AIF' button to draw an ROI, then click on 'Save ROI'.

4. ROI will be saved in MATLAB workspace as 'aif_roi'. Rename this variable to 'aifmask' and save it as 'AIF_Mask_P001GE_M.mat' in the folder that contains the 'P001' folder.

5. Draw an ROI for VOF and rename the new 'aif_roi' to 'veinmask'. And save it as 'Vein_Mask_P001GE_M.mat' in the same folder as above.


# Example notated Run 
Mira Liu April 26th 2023

Get data from the scanner, then run DicomSort.m on the downloaded DICOM directory. 
Then need to sort the perfusion code, and rename the folders. There is a python notebook I have for that (DSCPerfusionSorter_WIP.ipynb) but can also write it yourself or use Tim/Yong's code. 
Then... 

### Drawing AIF and Vein (VOF)

### Example Code
>> addpath '/Users/neuroimaging/Desktop/MR-Code/DSC_code/Perfusion/Code'\
>> addpath '/Users/neuroimaging/Desktop/MR-Code/DSC_code/Perfusion/Code/Other'\
>> addpath '/Users/neuroimaging/Desktop/MR-Code/DSC_code/Perfusion/Code/PerfCore'\
>> dcmpath = '/Users/neuroimaging/Desktop/DATA/Tabitha_Example/P001/';\
>> ep2d_dcmpath = [dcmpath, '/ep2d_perf'];\
>> dsc = loadDSC(ep2d_dcmpath);\
>> AIFviewer(dsc);

Now draw the AIF mask
>> AIFviewer(dsc);\
>> aifmask = aif_roi; %rename\
>> save ('/Users/neuroimaging/Desktop/AIF_Mask_P001GE_M.mat', 'aifmask') %savemask

Now also draw the Vein mask
>> AIFviewer(dsc);\
>> veinmask = aif_roi;\
>> save ('/Users/neuroimaging/Desktop/Vein_Mask_P001GE_M.mat', 'veinmask') %savemask

### Then Run qCBF post-processing
>> Auto_qCBF_Philips_DD_all('/Users/neuroimaging/Desktop/Data/Tabitha_Example/')

### White Matter Mask
Once the T1 map is done, exit and draw a white matter mask
>> load('/Users/neuroimaging/Desktop/DATA/Tabitha_Example/T1mapping/P001_T1map.mat')\
>> figure; imshow(images.t1.T1map_pre - images.t1.T1map_post,[0 50]),colormap(gca,'jet'),colorbar, truesize([500 500])\
>> WM_SS = roipoly;\
>> save ('/Users/neuroimaging/Desktop/WM_Mask_P001GE_M.mat', 'WM_SS')

### Then run final time
>> Auto_qCBF_Philips_DD_all('/Users/neuroimaging/Desktop/Data/Tabitha_Example/')

### Results Viewing
Results can be viewed with imshow, or if you'd like to view it as a volume can use imagestack. 
>> load(‘/Users/neuroimaging/Desktop/DATA/Tabitha_Example/Result_MSwcf2/P001GE_M.mat’)\
>> imagestack(images.DD.qCBF_SVD)
