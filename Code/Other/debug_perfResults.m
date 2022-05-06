function [ vargout ] = debug_perfResults( data, figs )
%debug_perfResults debug...
%   Detailed explanation goes here
%
% Author: Yong Ik Jeong
% Date: 2017-05-24

if nargin < 2
    figs = 0;
end
%cases = [3,4,5,6,7,9,10,11]';

totalperf = length(data);
totalcase = totalperf/2;
cases = zeros(totalcase,1);
count = 0;
divcase = 3;
for ii = 1:totalcase
    if mod(ii,divcase) == 1
        count = count + 1;
        aiffigs(count) = count;
        wmcbffigs(count) = count+divcase;
        wmcbffigs2(count) = count+divcase*2;
    end    
end

count = 0;
casecount = 0;
caseind = 0;
inj1ind = 0;
inj2ind = 0;
for ii = 1:totalperf
    id = data(ii).id;
    if mod(ii,2)
        casecount = casecount + 1;
        tmp = regexpi(id,'COLLAT_(\d\d)','tokens');
        cases(casecount) = str2double(tmp{1});
    end
    id = [id(1:strfind(id,'_')-1) '\' id(strfind(id,'_'):end)];
    date = data(ii).date;
    date = [date(1:strfind(date,'_')-1) '\' date(strfind(date,'_'):end)];
    perf = data(ii).perf;
    perf = [perf(1:strfind(perf,'_')-1) '\' perf(strfind(perf,'_'):end)];
    
    ROIs = data(ii).result.ROIs;
    cutoffs_ROIs = data(ii).result.cutoffs_ROIs;
    image_names = data(ii).result.image_names;
    images = data(ii).result.images;
    masks_ROIs = data(ii).result.masks_ROIs;
    
    if mod(ii,divcase*2) == 1
        count = count + 1;
    end
    
    % Plot AIFs
    if figs
        figure(aiffigs(count));
        subplot(divcase,2,mod(ii-1,6)+1);
        plot(ROIs.data.AIF);
        title([id '-' date '-' perf]);
    end
    
    bat = cutoffs_ROIs.AIF(2);
    [peak, peakx] = min(ROIs.data.AIF);
    hwidth = peakx-bat;
    slope = (peak-ROIs.data.AIF(bat))/hwidth;
    if mod(ii,2)
        inj1ind = inj1ind + 1;
        aifdata1(inj1ind,:) = [slope hwidth bat];
    else
        inj2ind = inj2ind + 1;
        aifdata2(inj2ind,:) = [slope hwidth bat];
    end      
    
    % Plot qCBF WM distr.
    slcnum = ROIs.positions.n_slice_WM_DSC;
    qcbf = images{strmatch('qCBF_nSVD',image_names)}(:,:,slcnum);
    wmmask = logical(masks_ROIs.WM_SS);
    bmask = logical(automaskns(qcbf));
    if mod(ii,2)
        wmmask2 = logical(data(ii+1).result.masks_ROIs.WM_SS);
        bmask2 = logical(automaskns(data(ii+1).result.images{strmatch('qCBF_nSVD',image_names)}(:,:,slcnum)));
    else
        wmmask2 = logical(data(ii-1).result.masks_ROIs.WM_SS);
        bmask2 = logical(automaskns(data(ii-1).result.images{strmatch('qCBF_nSVD',image_names)}(:,:,slcnum)));
    end
    jwmmask = logical(wmmask .* wmmask2);
    jbmask = logical(bmask .* bmask2);
    tmpmin = min(qcbf(jwmmask));
    tmpmax = max(qcbf(jwmmask));
    if figs
        figure(wmcbffigs(count));
        subplot(divcase,2,mod(ii-1,6)+1);
        hist(qcbf(jwmmask),tmpmin:(tmpmax-tmpmin)/20:tmpmax);
        title([id '-' date '-' perf]);
        %xlim([0 200]);
    end
    
    % Get qcbf mean
    %is tmpmin different from cutoff(1)?
    [ycbf,xcbf] = hist(qcbf(jwmmask),tmpmin:(tmpmax-tmpmin)/20:tmpmax);
    temp = find(ycbf);
    cutoff(1) = xcbf(min(temp));
    cutoff(2) = xcbf(max(temp));
%     tmpmin
%     cutoff(1)
    
    [a xmin] = min(abs(xcbf - cutoff(1)));
    [b xmax] = min(abs(xcbf - cutoff(2)));
    
    options = fitoptions('gauss1');
    options.Lower = [-Inf 0.01 -Inf];
    options.Upper = [Inf Inf Inf];
    
    fitteddata = fit(xcbf(xmin:xmax)',ycbf(xmin:xmax)','gauss1', options);
    if figs
        figure(wmcbffigs2(count));
        subplot(divcase,2,mod(ii-1,6)+1);
        plot(fitteddata, xcbf(xmin:xmax),ycbf(xmin:xmax))
        title([id '-' date '-' perf '-' num2str(fitteddata.b1)]);
        %xlim([0 200]);
    end

    %figure;imshow(jbmask);
    %figure;imshow(qcbf);
    if mod(ii,2)
        qcbfdata1(inj1ind,:) = fitteddata.b1;
        %yv(data(ii).result.images{10},[0 1200],'overlay',jwmmask)
        %qcbfdata1(inj1ind,:) = mean(qcbf(logical(jbmask)));
    else
        qcbfdata2(inj2ind,:) = fitteddata.b1;
        %qcbfdata2(inj2ind,:) = mean(qcbf(logical(jbmask)));
    end
    
    % T1 changes
    %dt1wm = 1000*( 1./mean(images{strmatch('T1map_post',image_names)}(jwmmask)) - 1./mean(images{strmatch('T1map_pre',image_names)}(jwmmask)) );
    %dt1ss = 1000*( 1./mean(images{strmatch('T1map_post',image_names)}(logical(masks_ROIs.SS))) - 1./mean(images{strmatch('T1map_pre',image_names)}(logical(masks_ROIs.SS))) );
    dt1wm = 1000*(1/ROIs.values.T1s.wm_post - 1/ROIs.values.T1s.wm_pre);
    dt1ss = 1000*(1/ROIs.values.T1s.blood_post - 1/ROIs.values.T1s.blood_pre);
    if mod(ii,2)
        %t1wmdata1(inj1ind,:) = [dt1wm mean(images{strmatch('T1map_pre',image_names)}(jwmmask)) mean(images{strmatch('T1map_post',image_names)}(jwmmask))];
        %t1ssdata1(inj1ind,:) = [dt1ss mean(images{strmatch('T1map_pre',image_names)}(logical(masks_ROIs.SS))) mean(images{strmatch('T1map_post',image_names)}(logical(masks_ROIs.SS)))];
        t1wmdata1(inj1ind,:) = [dt1wm ROIs.values.T1s.wm_pre ROIs.values.T1s.wm_post];
        t1ssdata1(inj1ind,:) = [dt1ss ROIs.values.T1s.blood_pre ROIs.values.T1s.blood_post];
    else
        %t1wmdata2(inj2ind,:) = [dt1wm mean(images{strmatch('T1map_pre',image_names)}(jwmmask)) mean(images{strmatch('T1map_post',image_names)}(jwmmask))];
        %t1ssdata2(inj2ind,:) = [dt1ss mean(images{strmatch('T1map_pre',image_names)}(logical(masks_ROIs.SS))) mean(images{strmatch('T1map_post',image_names)}(logical(masks_ROIs.SS)))];
        t1wmdata2(inj1ind,:) = [dt1wm ROIs.values.T1s.wm_pre ROIs.values.T1s.wm_post];
        t1ssdata2(inj1ind,:) = [dt1ss ROIs.values.T1s.blood_pre ROIs.values.T1s.blood_post];
    end
    
    if mod(ii,2)
        cbvsswm1(inj1ind,1) = ROIs.data.fitting.CBVssFast.b1;
    else
        cbvsswm2(inj2ind,1) = ROIs.data.fitting.CBVssFast.b1;
    end
    
    if mod(ii,2)
        cbvdscwm1(inj1ind,1) = ROIs.data.fitting.CBVdsc.b1;
    else
        cbvdscwm2(inj2ind,1) = ROIs.data.fitting.CBVdsc.b1;
    end
end

vargout = {aifdata1, aifdata2, t1wmdata1, t1wmdata2, t1ssdata1, t1ssdata2, qcbfdata1, qcbfdata2, cbvsswm1, cbvsswm2, cbvdscwm1, cbvdscwm2};

%fprintf('1slope \t 1hwidth \t 1bat \t 2slope \t 2hwdith \t 2bat\n');
fprintf('\tslope hwidth bat\n');
fprintf('%2d %9.4f %9.4f %9.4f  |  %9.4f %9.4f %9.4f\n',[cases aifdata1 aifdata2]');
fprintf('\n\tdR1wm mean_pre mean_post\n');
fprintf('%2d %9.4f %9.4f %9.4f  |  %9.4f %9.4f %9.4f\n',[cases t1wmdata1 t1wmdata2]');
fprintf('\n\tdR1ss mean_pre mean_post\n');
fprintf('%2d %9.4f %9.4f %9.4f  |  %9.4f %9.4f %9.4f\n',[cases t1ssdata1 t1ssdata2]');
fprintf('\n\tqcbf mean\n');
fprintf('%2d %9.4f  |  %9.4f\n',[cases qcbfdata1 qcbfdata2]');

end

