clc;clear all;close all;
%LOAD DATA
load('SD_AST_L1B.mat');%load infrared
load('SD_AST_L1B_T.mat');%load temperatureSD_AST_L1B = SD_AST_LIB;
SD_AST_L1B2 = SD_AST_L1B(11:end-10,110:670,:);
r1 = SD_AST_L1B2(:,:,2);
r2 = SD_AST_L1B2(:,:,3);
NDVI_Z = NDVI(r1,r2);
T2 = T(11:end-10,110:670,:);
% z1 = reshape(NDVI_Z,1,700*561);
% z2 = reshape(T(:,110:670,:),1,700*561);
% figure,plot(z1,z2,'*');
TVDI_Z = TVDI(NDVI_Z,T2);
DVI_Z = DVI(r1,r2);
% LAI_Z = LAI(r1,r2);
LAI_Z = LAI2(NDVI_Z);
RVI_Z = RVI(r1,r2);
VWC_Z = VWC(r1,r2);
SAVI_Z = SAVI(r1,r2);
for i = 1:14
    SD_AST_L1B3(:,:,i) = guiyihua(SD_AST_L1B2(:,:,i));
end
SWC_Z = SWC(SD_AST_L1B3,NDVI_Z);
figure,imagesc(NDVI_Z);
colormap(jet),title('NDVI'),colorbar;
figure,imagesc(DVI_Z);
colormap(jet),title('DVI'),colorbar;
figure,imagesc(LAI_Z);
colormap(jet),title('LAI'),colorbar;
figure,imagesc(RVI_Z);
colormap(jet),title('RVI'),colorbar;
figure,imagesc(VWC_Z);
colormap(jet),title('VWC'),colorbar;
figure,imagesc(SAVI_Z);
colormap(jet),title('SAVI'),colorbar;
figure,imagesc(SWC_Z);
colormap(jet),title('SWC'),colorbar;
figure,imagesc(TVDI_Z);
colormap(jet),title('TVDI'),colorbar;
coff_SD_L1B = SD_AST_L1B2(:,:,1:9);
coff_SD_L1B(:,:,10) = NDVI_Z;
coff_SD_L1B(:,:,11) = DVI_Z;
coff_SD_L1B(:,:,12) = LAI_Z;
coff_SD_L1B(:,:,13) = RVI_Z;
coff_SD_L1B(:,:,14) = VWC_Z;
coff_SD_L1B(:,:,15) = SAVI_Z;
coff_SD_L1B(:,:,16) = SWC_Z;
coff_SD_L1B(:,:,17) = TVDI_Z;


