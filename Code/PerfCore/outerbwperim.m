function pout = outerbwperim(image)

conn = conndef(2,'min');
mask = zeros(size(image));
temp = zeros(size(image));
mask(find(image))=1;

for i = 2:size(image,1)-1
    for j = 2:size(image,2)-1
        if image(i,j)
            temp(i-1,j-1)=1;temp(i,j-1)=1;temp(i+1,j-1)=1;
            temp(i-1,j  )=1;temp(i,j  )=1;temp(i+1,j  )=1;
            temp(i-1,j+1)=1;temp(i,j+1)=1;temp(i+1,j+1)=1;
        end
    end
end

            
if nargout==0
    figure;imshow(temp - mask)
else
    pout = temp - mask;
end            