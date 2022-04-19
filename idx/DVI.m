function DVI_Z = DVI(imagedata)
r1 = imagedata(:,:,3);
for i = 1:size(imagedata,1)
    for j = 1:size(imagedata,2)
        r2(i,j) = sum(imagedata(i,j,4:9))/6;
        DVI_Z(i,j) = r2(i,j)-r1(i,j);
    end
end
