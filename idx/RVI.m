function RVI_Z = RVI(imagedata)
r1 = imagedata(:,:,3);
for i = 1:size(imagedata,1)
    for j = 1:size(imagedata,2)
        r2(i,j) = sum(imagedata(i,j,4:9))/6;
        RVI_Z(i,j) = r2(i,j)/r1(i,j);
    end
end
RVI_Z(isnan(RVI_Z)) = 0;