function NDVI_Z = NDVI(imagedata)
r1 = imagedata(:,:,3);
for i = 1:size(imagedata,1)
    for j = 1:size(imagedata,2)
        r2(i,j) = sum(imagedata(i,j,4:9))/6;
        NDVI_Z(i,j) = (r2(i,j)-r1(i,j))/(r2(i,j)+r1(i,j));
    end
end
NDVI_Z(isnan(NDVI_Z)) = 0;
