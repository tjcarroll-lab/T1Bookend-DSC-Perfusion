function [ret] = make_3_by_3_BAT (images, image_names, patient_number,path,type)
%-------------------------------------------------------------------%
%                                                                   %
%                                                                   %
%-------------------------------------------------------------------%

raw_BAT  = images{strmatch('BATmap',image_names)};

[n,m]=size(raw_BAT(:,:,1));
[n,m,nslices]= size(raw_BAT);

total_pages = ceil(nslices/9);

%-------------------------------------------------------------------%
%  resize to 128 x 128                                              %
%-------------------------------------------------------------------%
BAT = zeros(128,128,nslices);
if ( n ~= 128)  
    for slice = 1:nslices
        rescale = 128/n;
        BAT(:,:,slice) = imresize(raw_BAT(:,:,slice),rescale,'bicubic');
    end
else
    BAT = raw_BAT;
end

if ( n < 128)  
    for slice = 1:nslices
        rescale = 128/n;
        BAT(:,:,slice) = imresize(raw_BAT(:,:,slice),rescale,'bicubic');
    end
end

%-----------------------------------------------------------------%
% mean                                                            %
%-----------------------------------------------------------------%
mean_BAT = 0.0;
n_voxels = 0.0;
for i=1:128;
for j=1:128
for k=1:nslices
   if (BAT(i,j,k) >0.001)
      mean_BAT = mean_BAT + BAT(i,j,k);
      n_voxels = n_voxels + 1;
   end
end
end
end
mean_BAT = mean_BAT/n_voxels
            
std_BAT = 0.0;
n_voxels = 0.0;
for i=1:128;
for j=1:128
for k=1:nslices
   if (BAT(i,j,k) >0.001)
   std_BAT = std_BAT + (BAT(i,j,k)-mean_BAT)*(BAT(i,j,k)-mean_BAT);
   n_voxels = n_voxels + 1;
   end
end
end
end
std_BAT = sqrt(std_BAT/n_voxels)

image_panel = zeros(384,384);

position.first = [1 129 257 ];
position.last  = [128 256 384];

beenby = 0;
for slice = 1:nslices
    
    loc_slice = mod( slice , 10);
    
    if ( slice > 9.5) 
        loc_slice = loc_slice + 1;
    end
   
    column = mod((loc_slice), 3);
    if( column == 0)
        column = 3;
    end
    row    = floor( (loc_slice-0.5)/3 )+1;

   page = floor(slice/9) + 1;

   image_panel( position.first(row):position.last(row), position.first(column):position.last(column)) = BAT(:,:,slice);

   %-----------------------------------------------------------------------
   % print page 1                                                         % 
   %-----------------------------------------------------------------------
   if ( page == 2 )
       if (beenby == 0 )
         h = figure(1);
         imshow(image_panel,[10 30]);colorbar, colormap jet; title( [ ' BAT sec ', patient_number, ...
                 ' Page ', int2str(page-1), ' of ',  int2str(total_pages)]);  
%          filename = [patient_number, ' Page ', int2str(page-1),'.tif'];
         filename = [path '\Maps\' patient_number '\' patient_number type '_BAT' int2str(page-1) '.tif'];
         text(100, 400, [' Mean ', num2str(mean_BAT), ' +/- ', num2str(std_BAT)]);
         saveas(h, filename);
         image_panel = zeros(384,384);
         beenby = 1;
       end
   end
   
end 

%-----------------------------------------------------------------------
% print page 2                                                         % 
%-----------------------------------------------------------------------
h = figure(2);
imshow(image_panel,[10 30]);colorbar, colormap jet; title( [ ' BAT sec ', patient_number,  ...
                 ' Page ', int2str(page-1), ' of ',  int2str(total_pages)]); 
% filename = [patient_number, ' Page ', int2str(page),'.tif'];
filename = [path '\Maps\' patient_number '\' patient_number type '_BAT' int2str(page) '.tif'];
saveas(h, filename);

return