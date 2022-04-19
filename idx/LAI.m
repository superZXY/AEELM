function LAI_Z = LAI(imagedata)
SAVI_Z = SAVI(imagedata);
LAI_Z = 585.19.*SAVI_Z.^4-582.*SAVI_Z.^2+204.8.*SAVI_Z-23.68;
