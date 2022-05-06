function resizeDICOM

path = pwd;
for lengthOfArray = 1:1500
    image = dicomread(strcat(num2str(lengthOfArray),'.dcm'));
    header = dicominfo(strcat(num2str(lengthOfArray),'.dcm'));
    header.Rows = 128;
    header.Columns = 128;
    image = imresize(image,128/112);
    dicomwrite(image,strcat(path,'/New/',num2str(lengthOfArray),'.dcm'),header);

end