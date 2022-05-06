function [CBVssFast_per, CBVssSlow_per,T1s] = CalcCBVSSmap(images,image_names,position,wmmask)

% predefine
t0 = 800;

% prepare
temp = images{strmatch('T1map_pre',image_names)};
mask = ones(size(temp));
mask(find(temp>=2000))=0;mask(find(temp<=0))=0;

T1map_pre   = fSpatialFilter(images{strmatch('T1map_pre',image_names,'exact')},[1 1],mask);
M0map_pre   = fSpatialFilter(images{strmatch('M0map_pre',image_names,'exact')},[1 1],mask);
InvFmap_pre = fSpatialFilter(images{strmatch('InvFmap_pre',image_names,'exact')},[1 1],mask);
T1map_post  = fSpatialFilter(images{strmatch('T1map_post',image_names,'exact')},[1 1],mask);
M0map_post  = fSpatialFilter(images{strmatch('M0map_post',image_names,'exact')},[1 1],mask);
InvFmap_post= fSpatialFilter(images{strmatch('InvFmap_post',image_names,'exact')},[1 1],mask);

%T1map_pre = savitzkyGolay2D_rle_coupling(224,224,T1map_pre,5,5,3);
%T1map_post = savitzkyGolay2D_rle_coupling(224,224,T1map_post,5,5,3);

%YIJ 20170615 smooth
% T1map_pre = imfilter(T1map_pre,fspecial('gaussian',[3 3],1));
% T1map_post = imfilter(T1map_post,fspecial('gaussian',[3 3],1));
M0map_pre(M0map_pre < 0) = 0;
M0map_post(M0map_post < 0) = 0;
% M0map_pre = imfilter(M0map_pre,fspecial('gaussian',[3 3],1));
% M0map_post = imfilter(M0map_post,fspecial('gaussian',[3 3],1));
% InvFmap_pre = imfilter(InvFmap_pre,fspecial('gaussian',[3 3],1));
% InvFmap_post = imfilter(InvFmap_post,fspecial('gaussian',[3 3],1));

CBVssFast_per = zeros(size(T1map_pre,1), size(T1map_pre,2));
CBVssSlow_per = zeros(size(T1map_pre,1), size(T1map_pre,2));

if isfield(position,'blood_pre') & isfield(position,'blood_post')
    T1_blood_pre = 0;    M0_blood_pre = 0;   InvF_blood_pre = 0;    tempcount = 0;
    T1_blood_post = 0;   M0_blood_post = 0;  InvF_blood_post = 0;
    for i = 1:length(position.blood_pre)
        if (T1map_pre(position.blood_pre{i}(1),position.blood_pre{i}(2)) ~= 0) & (T1map_pre(position.blood_pre{i}(1),position.blood_pre{i}(2)) < 3000) & isfinite(T1map_pre(position.blood_pre{i}(1),position.blood_pre{i}(2)))
            T1_blood_pre = T1_blood_pre + T1map_pre(position.blood_pre{i}(1),position.blood_pre{i}(2));
            tempcount = tempcount + 1;
        end
        M0_blood_pre = M0_blood_pre + M0map_pre(position.blood_pre{i}(1),position.blood_pre{i}(2));
        InvF_blood_pre = InvF_blood_pre + InvFmap_pre(position.blood_pre{i}(1),position.blood_pre{i}(2));
    end
    T1_blood_pre = T1_blood_pre/tempcount; 
    M0_blood_pre = M0_blood_pre/length(position.blood_pre);
    InvF_blood_pre = InvF_blood_pre/length(position.blood_pre);  
    
    tempcount = 0;
    for i = 1:length(position.blood_post)
        if (T1map_post(position.blood_post{i}(1),position.blood_post{i}(2)) ~= 0) & (T1map_post(position.blood_post{i}(1),position.blood_post{i}(2)) < 2000)
            T1_blood_post = T1_blood_post + T1map_post(position.blood_post{i}(1),position.blood_post{i}(2));
            tempcount = tempcount + 1;
        end
        M0_blood_post = M0_blood_post + M0map_post(position.blood_post{i}(1),position.blood_post{i}(2));
        InvF_blood_post = InvF_blood_post + InvFmap_post(position.blood_post{i}(1),position.blood_post{i}(2));
    end
    T1_blood_post = T1_blood_post/tempcount ; 
    M0_blood_post = M0_blood_post/length(position.blood_post) ; 
    InvF_blood_post = InvF_blood_post/length(position.blood_post) ; 
end

T1s.blood_pre = T1_blood_pre;
T1s.blood_post = T1_blood_post;
T1s.wm_pre = mean(T1map_pre(logical(wmmask)));
T1s.wm_post = mean(T1map_post(logical(wmmask)));
CBVssFast_per = 100*(1./T1map_pre - 1./T1map_post)/(1/T1_blood_pre - 1/T1_blood_post);

% norminator = M0map_pre.*(InvFmap_pre.*exp(-t0./T1map_pre)-1) - M0map_post.*(InvFmap_post.*exp(-t0./T1map_post)-1);
% denorminator = M0_blood_pre*(InvF_blood_pre*exp(-t0/T1_blood_pre)-1) - M0_blood_post*(InvF_blood_post*exp(-t0/T1_blood_post)-1);
norminator = (exp(-t0./T1map_pre)-1) - (exp(-t0./T1map_post)-1);
denorminator = (exp(-t0/T1_blood_pre)-1) - (exp(-t0/T1_blood_post)-1);
CBVssSlow_per = 100*norminator/denorminator;

% filtering
CBVssFast_per = fSpatialFilter(CBVssFast_per,[2 3],mask);
CBVssSlow_per = fSpatialFilter(CBVssSlow_per,[2 3],mask);



clear path_IR_pre image_names masks images position mask target count T t0 mask T1map_pre T1map_post M0map_pre M0map_post InvFmap_pre InvFmap_post T1_blood_pre T1_blood_post M0_blood_pre M0_blood_post InvF_blood_pre InvF_blood_post norminator denorminator

