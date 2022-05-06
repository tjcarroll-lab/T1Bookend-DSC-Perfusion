function fmask = AutoMaskNSforWM(images,threshold)

%%========================================================================%


mask = zeros(size(images,1),size(images,2));
mask(find(images(:,:,1)))=1;

%Exclude scalp

%% Method 1: Erosion and Dilation
mask1 = mask;
percentScalp = 0.80; 
brainArea = bwarea(mask);
%Reduce mask size 
index = 1;  count = 0;
while index
    temp = bwperim(mask);
    mask = mask - temp;
    maskarea = bwarea(mask);
    count = count + 1;
    if (maskarea < percentScalp*brainArea)
        index = 0;
    end
end
%Reincrease mask size 
while count
    temp = outerbwperim(mask);
    mask = mask + temp;
    count = count - 1;
end
%Final mask
fmask = mask.*mask1;


% %%%%%%Method 2: Ellipse drawing
% brainCover = bwarea(mask)/(size(images,1)*size(images,2));
% if (brainCover >= 0.25) && (brainCover < 0.30)
%     a = 0.55*size(images,1)/2;  %Major axis = 2a
%     b = 0.40*size(images,2)/2;  %Minor axis = 2b
% elseif (brainCover >= 0.30) && (brainCover < 0.35)
%     a = 0.60*size(images,1)/2;  
%     b = 0.45*size(images,2)/2;  
% elseif (brainCover >= 0.35) && (brainCover < 0.40)
%     a = 0.65*size(images,1)/2;  
%     b = 0.50*size(images,2)/2; 
% elseif (brainCover >= 0.40) && (brainCover < 0.45)
%     a = 0.70*size(images,1)/2;  
%     b = 0.55*size(images,2)/2;  
% else 
%     a = size(images,1);
%     b = size(images,2);
% end
% % % percentScalp = 0.20;
% % % center = [size(images,1)/2 size(images,2)/2];
% % % fmask_area = (1-percentScalp)*bwarea(mask);  
% % % ab = fmask_area/(4*pi);   %Ellipse area = 4*pi*a*b; b = 2/3*a;
% % % a = sqrt(3/2*ab);
% % % b = 2/3*a;
% fmask = mask;
% for i = 1:size(images,1)
%     for j = 1:size(images,2)
%         if ((i-center(1))^2/a^2 + (j-center(2))^2/b^2 > 1)     
%            fmask(i,j) = 0;     
%         end
%     end
% end


