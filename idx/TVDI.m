function TVDI = TVDI(NDVI_Z,T2)
% clc;clear all;close all;
% load('JS_AST_L1B.mat')
% load('JS_AST_L1B_T.mat')
% JS_AST_L1B2 = JS_AST_L1B(11:end-10,110:670,:);
% T2 = T(11:end-10,110:670,:);
% r1 = JS_AST_L1B2(:,:,2);
% r2 = JS_AST_L1B2(:,:,3);
% NDVI_Z = NDVI(r1,r2);

z1 = reshape(NDVI_Z,1,size(NDVI_Z,1)*size(NDVI_Z,2));
z2 = reshape(T2,1,size(T2,1)*size(T2,2));
figure,plot(z1,z2,'*');
cube_2(1,:) = z1(find(z1>=0&z1<0.2));
cube_2(2,:) = z2(find(z1>=0&z1<0.2));
% figure,plot(cube_2(1,:),cube_2(2,:),'*');
cube_4(1,:) = z1(find(z1>=0.2&z1<0.4));
cube_4(2,:) = z2(find(z1>=0.2&z1<0.4));
figure,plot(cube_4(1,:),cube_4(2,:),'*');
cube_6(1,:) = z1(find(z1>=0.4&z1<1));
cube_6(2,:) = z2(find(z1>=0.4&z1<1));
% figure,plot(cube_6(1,:),cube_6(2,:),'*');

% c1 = z1(find(z1>=0));
% c2 = z2(find(z1>=0));
c1 = cube_4(1,:);
c2 = cube_4(2,:);
[re_NDVI,pos_NDVI]=sort(c1,2,'ascend');%根据NDVI排序
T_pos = c2(pos_NDVI);
plot(re_NDVI,T_pos,'*');

a=T_pos;
b= re_NDVI;
% subplot(1,2,1)
% plot(a)
% axis tight

window=300;
window_3 = window/20;
num=floor(length(a)/window);
[in,out]=deal(zeros(1,num-1));
for i=1:num-1
    mysize=i*window+1:(i+1)*window;
    duan_T = a(mysize);
    duan_N = b(mysize);
    [re_duan_T,pos_duan_T]=sort(duan_T,2,'ascend');% 每一段根据温度重新排序
    re_duan_N = duan_N(pos_duan_T);
    in((i-1)*window_3+1:i*window_3)=re_duan_T(1:window_3);
    in_pos((i-1)*window_3+1:i*window_3)=re_duan_N(1:window_3);
    out((i-1)*window_3+1:i*window_3)=re_duan_T(end-window_3+1:end);
    out_pos((i-1)*window_3+1:i*window_3)=re_duan_N(end-window_3+1:end);
end
figure,plot(in_pos,in,'y*');
hold on;
plot(out_pos,out,'g*');
hold on;
f1 = fit(in_pos',in','poly1');        %低温度参数
plot(f1,in_pos,in);hold on;
f2 = fit(out_pos',out','poly1');      %高温度参数
plot(f2,out_pos,out);

TVDI = (T2-f1.p2-f1.p1*NDVI_Z)./((f2.p2+f2.p1*NDVI_Z)-(f1.p2-f1.p1*NDVI_Z));

% figure,imagesc(TDVI);
% colormap(jet),title('TDVI'),colorbar;
% tt1 = T2-f1.p2-f1.p1*NDVI_Z;
% tt2 = (f2.p2+f2.p1*NDVI_Z)-(f1.p2-f1.p1*NDVI_Z);
% figure,imagesc(tt1);
% colormap(jet),title('tt1'),colorbar;
% figure,imagesc(tt2);
% colormap(jet),title('tt2'),colorbar;
% figure,imagesc(NDVI_Z);
% colormap(jet),title('NDVI'),colorbar;
% figure,imagesc(T2);
% colormap(jet),title('T'),colorbar;