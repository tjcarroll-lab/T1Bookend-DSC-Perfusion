% CheckResults takes as an argument the patient folder number (eg, for P037
% you would call CheckResults(37)).
function CheckResults(patient_numbers, slice)
qCBFstorage = zeros(size(patient_numbers));
qCBF_WMarray = zeros(size(patient_numbers));
for patient = 1:length(patient_numbers)
    
    patient_number = patient_numbers(patient);
    if patient_number < 10
        patient_number = ['P00' num2str(patient_number)];
    elseif patient_number<100
        patient_number = ['P0' num2str(patient_number)];
    else
        patient_number = ['P' num2str(patient_number)];
    end
    
    %path = 'C:\Documents and Settings\CWloka\My Documents\MATLAB\MS Processing\BookEnd Processing\Main Directory';
    path = pwd;
    type = 'GE';
    load([path '\Result_MSwcf2\' patient_number type '_M.mat'])
    if(exist(strcat(path, '\Maps\', patient_number)) ~= 7)
        mkdir(strcat(path, '\Maps'), patient_number);
    end
    mean_qCBF = make_3_by_3_qCBF (images, image_names, patient_number,path,type)
    make_3_by_3_MTT (images, image_names, patient_number,path,type)
    make_3_by_3_qCBV (images, image_names, patient_number,path,type)
    qCBFstorage(patient) = mean_qCBF;
    
    h = figure(10);
    preT1 = images{strmatch('T1map_pre', image_names)};
    imshow(preT1, [0 1200]); colormap jet; colorbar;
    filename = [path '\Maps\' patient_number '\' patient_number type '_preT1.tif'];
    saveas(h, filename); close(h);
    h = figure(10);
    postT1 = images{strmatch('T1map_post', image_names)};
    imshow(postT1, [0 1200]); colormap jet; colorbar;
    filename = [path '\Maps\' patient_number '\' patient_number type '_postT1.tif'];
    saveas(h, filename); close(h);
    h = figure(10);
    preMinusPostT1 = preT1-postT1;
    imshow(preMinusPostT1, [-25 125]); colormap jet; colorbar;
    filename = [path '\Maps\' patient_number '\' patient_number type '_preMinusPostT1.tif'];
    saveas(h, filename); close(h);
    
    qCBF = images{strmatch('qCBF_nSVD',image_names)};
    [WM_qCBF std] = fast_roi(masks_ROIs.WM_SS,qCBF(:,:,ROIs.positions.n_slice_WM_DSC))
    qCBF_WMarray(patient) = WM_qCBF;
    qCBV = images{strmatch('qCBV_DSC',image_names)};
    [WM_qCBV std] = fast_roi(masks_ROIs.WM_SS,qCBV(:,:,ROIs.positions.n_slice_WM_DSC))

    
    %GM qCBF
    figure
    imshow(qCBF(:,:,slice),[0 120]);colormap jet; colorbar;
    %GMroi = roipoly;
    %save(['GMrois\' patient_number '_GMroi.mat'],'GMroi','slice')
    % load(['GMrois\' patient_number '_GMroi.mat'],'GMroi','slice')
    %[GM_qCBF std] = fast_roi(GMroi,qCBF(:,:,slice))
    
%     
%     figure(3)
%     imshow(images{strmatch('T1map_pre',image_names)},[0 1000]);colormap jet; colorbar;
%     title('T1 map pre')
%     figure(4)
%     imshow(images{strmatch('T1map_post',image_names)},[0 1000]);colormap jet; colorbar;
%     title('T1 map post')
    
    ROIs.values.T1s
    
    T1map_pre = images{strmatch('T1map_pre',image_names)};
    T1map_post = images{strmatch('T1map_post',image_names)};
    [WMmask_T1pre std]= fast_roi(masks_ROIs.WM_SS,T1map_pre)
    [WMmaks_T1post std]= fast_roi(masks_ROIs.WM_SS,T1map_post)
    
    ROIs.data.fitting.CBVdsc
    ROIs.data.fitting.CBVssFast
    ROIs.values.CF
%     
%     figure(5)
%     imshow(masks_ROIs.SS)
%     title('blood voxels selected')
%     figure(6)
%     imshow(masks_ROIs.WM_SS)
%     title('WM mask')
%     figure(7)
%     plot(ROIs.data.AIF)
%     title('AIF signal')
end
save('qCBFstorage');
end





