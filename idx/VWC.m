function VWC_Z = VWC(imagedata)
NDVI1 = NDVI(imagedata);
VWC_Z = 192.64.*NDVI1.^5-417.46.*NDVI1.^4+347.96.*NDVI1.^3-138.00.*NDVI1.^2+30.699.*NDVI1-2.82;
