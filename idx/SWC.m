function SWC_Z = SWC(imagedata)
A1 = imagedata(:,:,1);
A2 = imagedata(:,:,2);
A3 = imagedata(:,:,3);
SWC_Z = 0.77-42.19.*log10((0.6968.*A1+0.5228.*A2-0.2237.*A3+20.20)./...
    (1.089-0.00579.*A3+0.003308.*A1+0.002482.*A2)-18.0);
SWC_Z(isnan(SWC_Z)) = 0;