%% I. 清空环境变量
clear all
clc
close all
%% II. 训练集/测试集产生
%%
% 1. 导入数据
% load train_g22.mat;
load coff_SD_L1B.mat;%load infrared data
load('SD_AST_L1B_T.mat');%load land surface temperature
[coff_x, coff_y, coff_z] = size(coff_SD_L1B);
reshape_coff = reshape(coff_SD_L1B,coff_x*coff_y, coff_z);
reshape_T = reshape(T(11:end-10,110:670,:),coff_x*coff_y,1);
figure,imagesc(T(11:end-10,110:670,:));
colormap(jet),title('温度反演结果'),colorbar;
clear T;
pos5 = reshape_coff';
tol5 = reshape_T';
% 测试集——个样本


NDVI_line = pos5(10,:);
NDVI_f_0_2 = find(NDVI_line<-0.2);
NDVI_0 = find(NDVI_line>=-0.2 & NDVI_line<0);
NDVI_0_2 = find(NDVI_line>=0 & NDVI_line<0.2);
% NDVI_0_5 = find(NDVI_line>=0.2);

NDVI_0_5 = find(NDVI_line>=0.2 );
% NDVI_1 = find(NDVI_line>=0.5);

T_test = tol5;
P_test1 = pos5([1:10,15,17],NDVI_f_0_2);
T_test1 = tol5(NDVI_f_0_2);
P_test2 = pos5([1:11,13,16,17],NDVI_0);
T_test2 = tol5(NDVI_0);
P_test3 = pos5([1:9,16,17],NDVI_0_2);
T_test3 = tol5(NDVI_0_2);
P_test4 = pos5([1:10,12,16,17],NDVI_0_5);
T_test4 = tol5(NDVI_0_5);
% P_test5 = pos5([1:10,12,13,17],NDVI_1);
% T_test5 = tol5(NDVI_1);
%%
% nc = 11;
k=1;gate1=1;gate2=1;gate3=1;gate4=1;gate5=1;
% 2. 随机产生训练集和测试集
% for ii = 1:nc
ii=1;
while 1
    
    temp1 = randperm(size(pos5(NDVI_f_0_2),2));
    temp2 = randperm(size(pos5(NDVI_0),2));
    temp3 = randperm(size(pos5(NDVI_0_2),2));
    temp4 = randperm(size(pos5(NDVI_0_5),2));
%     temp5 = randperm(size(pos5(NDVI_1),2));
    % 训练集——1200个样本
    %%
    P_train1 = pos5([1:10,15,17],temp1(1:sqrt(size(temp1,2))));%NDVI<-0.2
    T_train1 = tol5(:,temp1(1:sqrt(size(temp1,2))));
    P_train2 = pos5([1:11,13,16,17],temp2(1:sqrt(size(temp2,2))));%-0.2<NDVI<0
    T_train2 = tol5(:,temp2(1:sqrt(size(temp2,2))));
    P_train3 = pos5([1:9,16,17],temp3(1:sqrt(size(temp3,2))));%0<NDVI<0.2
    T_train3 = tol5(:,temp3(1:sqrt(size(temp3,2))));
    P_train4 = pos5([1:10,12,16,17],temp4(1:sqrt(size(temp4,2))));%0.2<NDVI<0.5
    T_train4 = tol5(:,temp4(1:sqrt(size(temp4,2))));
%     P_train5 = pos5([1:10,12,13,17],temp5(1:sqrt(size(temp5,2))));%0.5<NDVI<1
%     T_train5 = tol5(:,temp5(1:sqrt(size(temp5,2))));
    %% III. 数据归一化
    %%
    % 1. 训练集
    [Pn_train1,inputps1] = mapminmax(P_train1);
    Pn_test1 = mapminmax('apply',P_test1,inputps1);
    [Pn_train2,inputps2] = mapminmax(P_train2);
    Pn_test2 = mapminmax('apply',P_test2,inputps2);
    [Pn_train3,inputps3] = mapminmax(P_train3);
    Pn_test3 = mapminmax('apply',P_test3,inputps3);
    [Pn_train4,inputps4] = mapminmax(P_train4);
    Pn_test4 = mapminmax('apply',P_test4,inputps4);
%     [Pn_train5,inputps5] = mapminmax(P_train5);
%     Pn_test5 = mapminmax('apply',P_test5,inputps5);
    %%
    % 2. 测试集
    [Tn_train1,outputps1] = mapminmax(T_train1);
    Tn_test1 = mapminmax('apply',T_test1,outputps1);
    [Tn_train2,outputps2] = mapminmax(T_train2);
    Tn_test2 = mapminmax('apply',T_test2,outputps2);
    [Tn_train3,outputps3] = mapminmax(T_train3);
    Tn_test3 = mapminmax('apply',T_test3,outputps3);
    [Tn_train4,outputps4] = mapminmax(T_train4);
    Tn_test4 = mapminmax('apply',T_test4,outputps4);
%     [Tn_train5,outputps5] = mapminmax(T_train5);
%     Tn_test5 = mapminmax('apply',T_test5,outputps5);

    %% IV. ELM创建/训练
    [Beta01,out_w1,IW1,B1,LW1,TF1,TYPE1] = elmtrain(Pn_train1,Tn_train1,floor(log2(size(Pn_train1,2))),'softplus',0);
    [Beta02,out_w2,IW2,B2,LW2,TF2,TYPE2] = elmtrain(Pn_train2,Tn_train2,floor(log2(size(Pn_train2,2))),'softplus',0);
    [Beta03,out_w3,IW3,B3,LW3,TF3,TYPE3] = elmtrain(Pn_train3,Tn_train3,floor(log2(size(Pn_train3,2))),'softplus',0);
    [Beta04,out_w4,IW4,B4,LW4,TF4,TYPE4] = elmtrain(Pn_train4,Tn_train4,floor(log2(size(Pn_train4,2))),'softplus',0);
%     [Beta05,out_w5,IW5,B5,LW5,TF5,TYPE5] = elmtrain(Pn_train5,Tn_train5,floor(log2(size(Pn_train5,2))),'softplus',0);
    %% V. ELM仿真测试
    tn_sim1 = elmpredict(Pn_test1,IW1,B1,LW1,TF1,TYPE1);
    tn_sim2 = elmpredict(Pn_test2,IW2,B2,LW2,TF2,TYPE2);
    tn_sim3 = elmpredict(Pn_test3,IW3,B3,LW3,TF3,TYPE3);
    tn_sim4 = elmpredict(Pn_test4,IW4,B4,LW4,TF4,TYPE4);
%     tn_sim5 = elmpredict(Pn_test5,IW5,B5,LW5,TF5,TYPE5);

    %%
    % 1. 反归一化
    T_sim1(ii,:) = mapminmax('reverse',tn_sim1,outputps1);
    T_sim2(ii,:) = mapminmax('reverse',tn_sim2,outputps2);
    T_sim3(ii,:) = mapminmax('reverse',tn_sim3,outputps3);
    T_sim4(ii,:) = mapminmax('reverse',tn_sim4,outputps4);
%     T_sim5(ii,:) = mapminmax('reverse',tn_sim5,outputps5);
    
    if ii>1 & mod(ii,2) ~= 0
        if gate1 == 1
        [rank_pix1,I1] = sort(T_sim1);
        result_T1(k,:) = rank_pix1((ii-1)/2+1,:);
        end
        if gate2 == 1
        [rank_pix2,I2] = sort(T_sim2);
        result_T2(k,:) = rank_pix2((ii-1)/2+1,:);
        end
        if gate3 == 1
        [rank_pix3,I3] = sort(T_sim3);
        result_T3(k,:) = rank_pix3((ii-1)/2+1,:);
        end
        if gate4 == 1
        [rank_pix4,I4] = sort(T_sim4);
        result_T4(k,:) = rank_pix4((ii-1)/2+1,:);
        end
%         if gate5 == 1
%         [rank_pix5,I1] = sort(T_sim5);
%         result_T5(k,:) = rank_pix5((ii-1)/2+1,:);
%         end
        
        if k>2
            if gate1 == 1
                RMSE1 = sqrt(mse(result_T1(k,:) - result_T1(k-1,:)));
            end
            if gate2 == 1
                RMSE2 = sqrt(mse(result_T2(k,:) - result_T2(k-1,:)));
            end
            if gate3 == 1
                RMSE3 = sqrt(mse(result_T3(k,:) - result_T3(k-1,:)));
            end
            if gate4 == 1
                RMSE4 = sqrt(mse(result_T4(k,:) - result_T4(k-1,:)));
            end
%             if gate5 == 1
%                 RMSE5 = sqrt(mse(result_T5(k,:) - result_T5(k-1,:)));
%             end
            disp(['RMSE1',num2str(ii),' = ',num2str(RMSE1)]);
            disp(['RMSE2',num2str(ii),' = ',num2str(RMSE2)]);
            disp(['RMSE3',num2str(ii),' = ',num2str(RMSE3)]);
            disp(['RMSE4',num2str(ii),' = ',num2str(RMSE4)]);
%             disp(['RMSE5',num2str(ii),' = ',num2str(RMSE5)]);
            if RMSE1 < 0.01
                gate1 = 2;
            end
            if RMSE2 < 0.01
                gate2 = 2;
            end
            if RMSE3 < 0.01
                gate3 = 2;
            end
            if RMSE4 < 0.01
                gate4 = 2;
            end
%             if RMSE5 < 0.01
%                 gate5 = 2;
%             end
            if gate1*gate2*gate3*gate4 ==16
                break;
            end
        end
        k = k+1;
       
    end
     ii=ii+1;
    clear Tn_test Tn_train Pn_train Pn_test P_train T_train ;
end


T_sim(:,NDVI_f_0_2) = result_T1(end,:);
T_sim(:,NDVI_0) = result_T2(end,:);
T_sim(:,NDVI_0_2) = result_T3(end,:);
T_sim(:,NDVI_0_5) = result_T4(end,:);
% T_sim(:,NDVI_1) = result_T5(end,:);
T_sim(T_sim<0)=0;
% for ii = 1:nc
%     T(:,:,ii) = reshape(T_sim(ii,:),coff_x,coff_y);
% end
result_T = reshape(T_sim,coff_x,coff_y);
% for i = 1:coff_x
%     for j = 1:coff_y
%         pix = reshape(T(i,j,:),1,nc);
%         [rank_pix,I] = sort(pix);
%         result_T(i,j) = rank_pix((nc-1)/2+1);
% %         result_T(i,j) = rank_pix(nc);
%         clear pix rank_pix I;
%     end
% end
for i = 1:size(result_T,1)
    for j = 1:size(result_T,2)
        if result_T(i,j) < 0 || result_T(i,j) > 42
            result_T(i,j) = 0;
        end
    end
end
figure,imagesc(result_T);
colormap(jet),title('温度预测结果'),colorbar;
% save T_sim.mat T_sim
%% VI. 结果对比
clear T_sim;
T_sim = reshape(result_T,1,coff_x*coff_y);
result = [T_test' T_sim'];

%%
% 1. 相对误差error
error = abs(T_sim - T_test)./T_test;
error2 = abs(T_sim - T_test);



%%
% 1. 均方误差
E = mse(T_sim - T_test);
MSE = sqrt(mse(T_sim - T_test));
%%
% 2. 决定系数
N = length(T_test);
R2=(N*sum(T_sim.*T_test)-sum(T_sim)*sum(T_test))^2/((N*sum((T_sim).^2)-(sum(T_sim))^2)*(N*sum((T_test).^2)-(sum(T_test))^2)); 

%% VII. 绘图
figure,plot(error);
xlabel('x axis');
ylabel('y axis');
legend('error');
xlabel('预测样本')
ylabel('误差')
string = {'力矩预测误差(ELM)'};
title(string)

figure,plot(1:N,T_test,'r-*',1:N,T_sim,'b:o')
% T_sim(find(T_sim>40))=0;
% plot(T_test',T_sim','b:o')

grid on
legend('真实值','预测值')
xlabel('预测样本')
ylabel('力矩值')
string = {'测试集力矩预测结果对比(ELM)';['(mse = ' num2str(E) ' R^2 = ' num2str(R2) ')']};
title(string)


% imshow(aa(:,:,1),[]);
% aa = reshape(aa,1,size(aa,1)*size(aa,2));
aa = reshape(T_sim,coff_x,coff_y);
% aa = reshape(T_sim,300,310);
for i = 1:size(aa,1)
    for j = 1:size(aa,2)
        if aa(i,j) < 2 || aa(i,j) > 42
            aa(i,j) = 0;
        end
    end
end

cut_pre = aa;  %500*500预测温度方阵
cut_pre_line = reshape(cut_pre,size(cut_pre,1)*size(cut_pre,2),1);   %温度预测值
true = reshape(tol5,coff_x,coff_y);
sqr_true = true;
true_line =  reshape(sqr_true,size(sqr_true,1)*size(sqr_true,2),1); %温度真实值
figure,plot(cut_pre_line,true_line,'*');
xlabel('预测值')
ylabel('真实值')
figure,imshow(sqr_true,[]);
figure,imshow(cut_pre,[]);
cc = abs(cut_pre-sqr_true);
figure,imagesc(cc);
colormap(jet),colorbar;
