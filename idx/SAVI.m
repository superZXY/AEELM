function SAVI_Z = SAVI(imagedata)
L = 0.5;
r1 = imagedata(:,:,3);
for i = 1:size(imagedata,1)
    for j = 1:size(imagedata,2)
        r2(i,j) = sum(imagedata(i,j,4:9))/6;
        SAVI_Z(i,j) = (1+L)*(r2(i,j)-r1(i,j))/(r2(i,j)+r1(i,j)+L);
    end
end
SAVI_Z(isnan(SAVI_Z)) = 0;