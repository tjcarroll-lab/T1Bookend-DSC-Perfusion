function [position_SS,SSmask,info] = AutoSS(varargin)
%function [position_SS,perCutoff] = AutoSS(T1map_pre,T1map_post,M0map_pre,M0map_post)

% Automation of Saggittal Sinus detection 
% modified (WYS) 12/19/2006

%% define varaibles
T1map_pre  = varargin{1};
T1map_post = varargin{2};
M0map_pre  = varargin{3};
M0map_post = varargin{4};
RESNORMmap_pre = varargin{5};
RESNORMmap_post = varargin{6};
info = cell(0);
global injectionNum;
global glblTargetPath;
global seqType;

%%%%Grady adjustment to simply draw SS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dR1map = 1000*(1./T1map_post - 1./T1map_pre);  % sec
dR1map(find(~isfinite(dR1map)))=0;

%% make a mask
%mask = AutoMaskNSforWM(T1map_pre);
%YIJ 20170614 smooth mask
mask = automaskns(imfilter(T1map_pre,fspecial('gaussian',[3 3],2)));
mask = imfill(mask);
% T1map_pre = imfilter(T1map_pre,fspecial('gaussian',[3 3],1));
% T1map_post = imfilter(T1map_post,fspecial('gaussian',[3 3],1));
M0map_pre(M0map_pre < 0) = 0;
M0map_post(M0map_post < 0) = 0;
% M0map_pre = imfilter(M0map_pre,fspecial('gaussian',[3 3],1));
% M0map_post = imfilter(M0map_post,fspecial('gaussian',[3 3],1));

% YIJ 20170509: draw mask of dog brain
%f1 = figure; imshow(T1map_pre, []); title('Make mask for autoSS'); mask = roipoly;close(f1);

mask(find(T1map_pre<=0)) =0;
mask(find(T1map_post<=0)) =0;
%YIJ 20170608
mask(find(T1map_pre>=2000)) = 0;
mask(find(T1map_post>=2000)) = 0;

%YIJ 20171013
mask(find(RESNORMmap_pre >= 10000)) = 0;
mask(find(RESNORMmap_post >= 10000)) = 0;

MeanT1pre  = mean(T1map_pre(find(mask)));
StdsT1pre  = std(T1map_pre(find(mask)));
MeanT1post = mean(T1map_post(find(mask)));
StdsT1post = std(T1map_post(find(mask)));

% make a baseline using M0
sumM0map = (M0map_pre + M0map_post).*mask;
sumM0map(sumM0map < 0) = 0;
MaxM0 = max(max(sumM0map));

[histM0,  M0table]  = hist(sumM0map(find(mask)),0:1:MaxM0);
tmp = histM0 == 0;
histM0(tmp) = [];
M0table(tmp) = [];
cutoffs_M0 = length(M0table);

loopingindex = 1;
preventEndless = 0;
while loopingindex
    tmask = mask;
    %tmask(1:size(mask,1)/2,:)=0;
    tmask(find(sumM0map<=M0table(cutoffs_M0)))=0;
    
    % check
    tmask(find(T1map_pre<MeanT1pre+StdsT1pre))=0;
    tmask(find(T1map_post>MeanT1post-StdsT1post))=0;

    count = sum(sum(tmask));
    %YIJ 20170614 include more points
    if (count > 10) | (cutoffs_M0 < 2)
        SSmask = tmask;
        loopingindex = 0;
    end
    cutoffs_M0 = cutoffs_M0 -1;
    preventEndless = preventEndless + 1;
    if preventEndless > 1000;
        loopingindex = 0;
        SSmask = tmask;
    end
end

%YIJ 20170614 higher threshold
perCutoff_T1pre  = 10; % 0.1
perCutoff_T1post = 10; % 0.1

%YIJ 20200723 remove this for now, because matrix size can be different
%between injections
%YIJ 20170614 include 1st inj ss mask 
% if injectionNum > 1
%     disp(['injection num > 1']);
%     if exist([glblTargetPath '/AutoSS/P001_AutoSS.mat'],'file')
%         matobj = matfile([glblTargetPath '/AutoSS/P001_AutoSS.mat']);
%         tempMask = matobj.('masks_ROIs');
%         tmpSSmask = SSmask;
%         tmpSSmask(logical(tempMask.SS)) = 1;
%         
%         meanT1pre  = mean(T1map_pre(find(SSmask)));
%         stdsT1pre  = std(T1map_pre(find(SSmask)));
%         meanT1post = mean(T1map_post(find(SSmask)));
%         stdsT1post = std(T1map_post(find(SSmask)));
%         nstdpre1 = stdsT1pre/meanT1pre;
%         nstdpost1 = stdsT1post/meanT1post;
%         
%         tmpmeanT1pre  = mean(T1map_pre(find(tmpSSmask)));
%         tmpstdsT1pre  = std(T1map_pre(find(tmpSSmask)));
%         tmpmeanT1post = mean(T1map_post(find(tmpSSmask)));
%         tmpstdsT1post = std(T1map_post(find(tmpSSmask)));
%         tmpnstdpre1 = tmpstdsT1pre/tmpmeanT1pre;
%         tmpnstdpost1 = tmpstdsT1post/tmpmeanT1post;
%         
%         %if tmpnstdpre1 < nstdpre1 && tmpnstdpost1 < nstdpost1
%         SSmask = tmpSSmask;
%         %end
%     end
% end

meanT1pre  = mean(T1map_pre(find(SSmask)));
stdsT1pre  = std(T1map_pre(find(SSmask)));
meanT1post = mean(T1map_post(find(SSmask)));
stdsT1post = std(T1map_post(find(SSmask)));
nstdpre1 = stdsT1pre/meanT1pre;
nstdpost1 = stdsT1post/meanT1post;

%YIJ 20170614 this method seems to work better
tmask = SSmask;
%tmask1 = SSmask;
%tmask2 = SSmask;
preventEndless = 0;
looppre = 1;
looppost = 1;
while looppre || looppost
    if 100*nstdpre1 > perCutoff_T1pre
        tmask1 = tmask;
        tmask2 = tmask;
        tmask1(find(T1map_pre>meanT1pre+stdsT1pre))=0;
        tmask2(find(T1map_pre<meanT1pre-stdsT1pre))=0;       
        tmeanT1pre1  = mean(T1map_pre(find(tmask1)));
        tstdsT1pre1  = std(T1map_pre(find(tmask1)));
        tnstdpre21 = tstdsT1pre1/tmeanT1pre1;
        tmeanT1pre2  = mean(T1map_pre(find(tmask2)));
        tstdsT1pre2  = std(T1map_pre(find(tmask2)));
        tnstdpre22 = tstdsT1pre2/tmeanT1pre2;
        if tnstdpre21 < tnstdpre22
            tmask = tmask1;
        else
            tmask = tmask2;
        end
    else
        looppre = 0;
    end
    if 100*nstdpost1 > perCutoff_T1post
        tmask1 = tmask;
        tmask2 = tmask;
        tmask1(find(T1map_post>meanT1post+stdsT1post))=0;
        tmask2(find(T1map_post<meanT1post-stdsT1post))=0;       
        tmeanT1post1  = mean(T1map_post(find(tmask1)));
        tstdsT1post1  = std(T1map_post(find(tmask1)));
        tnstdpost21 = tstdsT1post1/tmeanT1post1;
        tmeanT1post2  = mean(T1map_post(find(tmask2)));
        tstdsT1post2  = std(T1map_post(find(tmask2)));
        tnstdpost22 = tstdsT1post2/tmeanT1post2;
        if tnstdpost21 < tnstdpost22
            tmask = tmask1;
        else
            tmask = tmask2;
        end
    else       
        looppost = 0;
    end
    
    meanT1pre  = mean(T1map_pre(find(tmask)));
    stdsT1pre  = std(T1map_pre(find(tmask)));
    nstdpre2 = stdsT1pre/meanT1pre;
    meanT1post  = mean(T1map_post(find(tmask)));
    stdsT1post  = std(T1map_post(find(tmask)));
    nstdpost2 = stdsT1post/meanT1post;
    
    if nstdpre1 == nstdpre2 && nstdpost1 == nstdpost2
        break
    else
        nstdpre1 = nstdpre2;
        nstdpost1 = nstdpost2;
    end
    preventEndless = preventEndless + 1;
    if preventEndless > 100;
        nstdpre1 = 0;
        nstdpost1 = 0;
    end
end 

%YIJ 20170614 Original method commented out
% tmask1 = SSmask;
% preventEndless = 0;
% while (100*nstdpre1 > perCutoff_T1pre) 
%     tmask1(find(T1map_pre>meanT1pre+stdsT1pre))=0;
%     tmask1(find(T1map_pre<meanT1pre-stdsT1pre))=0;
%     meanT1pre  = mean(T1map_pre(find(tmask1)));
%     stdsT1pre  = std(T1map_pre(find(tmask1)));
%     nstdpre2 = stdsT1pre/meanT1pre;
%     if nstdpre1 == nstdpre2
%         break
%     else
%         nstdpre1 = nstdpre2;
%     end
%     preventEndless = preventEndless + 1;
%     if preventEndless > 100;
%         nstdpre1 = 0;
%     end
% end
% 
% tmask2 = SSmask;
% preventEndless = 0;
% while (100*nstdpost1 > perCutoff_T1post) 
%     tmask2(find(T1map_post>meanT1post+stdsT1post))=0;
%     tmask2(find(T1map_post<meanT1post-stdsT1post))=0;
%     meanT1post  = mean(T1map_post(find(tmask2)));
%     stdsT1post  = std(T1map_post(find(tmask2)));
%     nstdpost2 = stdsT1post/meanT1post;
%     if nstdpost1 == nstdpost2
%         break
%     else
%         nstdpost1 = nstdpost2;
%     end
%     preventEndless = preventEndless + 1;
%     if preventEndless > 100;
%         nstdpost1 = 0;
%     end
% end    

%final
%tmask = tmask1.*tmask2;

SSmask = tmask;

%if 2nd or higher injection use combined SS ROIs;
% if injectionNum > 1
%     disp(['injection num > 1']);
%     matobj = matfile([glblTargetPath '/AutoSS/P001_AutoSS.mat']);
%     tempMask = matobj.('masks_ROIs');
%     SSmask(logical(tempMask.SS)) = 1;    
% end
%fprintf('blood_pre: %f, blood_post: %f\n',mean(T1map_pre(logical(SSmask))),mean(T1map_post(logical(SSmask))));

% must be more than 2 pixels
if (sum(sum(SSmask)) < 2) 
    % 20200707 YIJ: Use vein mask instead
    if exist([glblTargetPath '\Vein_Mask_P' sprintf('%03d',injectionNum) seqType '.mat'],'file')
        tmpload = load([glblTargetPath '\Vein_Mask_P' sprintf('%03d',injectionNum) seqType '.mat']);
        SSmask = tmpload.veinmask;
        warning('less than 2 pixel SS');
        fprintf('AutoSS failed. Using vein mask instead\n');
    else
        error('less than 2 pixel SS');
    end
    %SSmask = T1map_pre > 1000 & T1map_pre < 1900 & T1map_post > 200 & T1map_post < 700;
    %disp('AutoSS failed. Made SSmask based on hard threshold');
    %f1 = figure; imshow(dR1map, [0 3]); colormap jet; title('(autoSS failed) Make roi of SS'); SSmask = roipoly;close(f1);
end

T1map_pre(logical(SSmask))
T1map_post(logical(SSmask))
fprintf('blood_pre: %f, blood_post: %f\n',mean(T1map_pre(logical(SSmask))),mean(T1map_post(logical(SSmask))));

%%%End grady Adjustment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Kevin modification of Grady adjustment 10/14/2016, did delta T1map in order to more easily locate
%the SS
%     figure;imshow(T1map_post-T1map_pre, []);
% %
%     %imshow(T1map_pre,[]);
%     SSmask = roipoly;
% 

position_SS = cell(0);
% save data
[X Y] = find(SSmask);
for i = 1:length(X)
	position_SS{i} = [X(i) Y(i)];
end
