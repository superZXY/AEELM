function [Beta0,out_w,IW,B,LW,TF,TYPE] = elmtrain(P,T,N,TF,TYPE)
% ELMTRAIN Create and Train a Extreme Learning Machine
% Syntax
% [IW,B,LW,TF,TYPE] = elmtrain(P,T,N,TF,TYPE)
% Description
% Input
% P   - Input Matrix of Training Set  (R*Q)
% T   - Output Matrix of Training Set (S*Q)
% N   - Number of Hidden Neurons (default = Q)
% TF  - Transfer Function:
%       'sig' for Sigmoidal function (default)
%       'sin' for Sine function
%       'hardlim' for Hardlim function
% TYPE - Regression (0,default) or Classification (1)
% Output
% IW  - Input Weight Matrix (N*R)
% B   - Bias Matrix  (N*1)
% LW  - Layer Weight Matrix (N*S)
% Example£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨£¨
% Regression:
% [IW,B,LW,TF,TYPE] = elmtrain(P,T,20,'sig',0)
% Y = elmtrain(P,IW,B,LW,TF,TYPE)
% Classification
% [IW,B,LW,TF,TYPE] = elmtrain(P,T,20,'sig',1)
% Y = elmtrain(P,IW,B,LW,TF,TYPE)
% See also ELMPREDICT
% Yu Lei,11-7-2010
% Copyright www.matlabsky.com
% $Revision:1.0 $
if nargin < 2
    error('ELM:Arguments','Not enough input arguments.');
end
if nargin < 3
    N = size(P,2);
end
if nargin < 4
    TF = 'sig';
end
if nargin < 5
    TYPE = 0;
end
if size(P,2) ~= size(T,2)
    error('ELM:Arguments','The columns of P and T must be same.');
end
[R,Q] = size(P);
if TYPE  == 1
    T  = ind2vec(T);
end
[S,Q] = size(T);
% Randomly Generate the Input Weight Matrix
IW = rand(N,R) * 2 - 1;
% Randomly Generate the Bias Matrix
B = rand(N,1);
BiasMatrix = repmat(B,1,Q);
% Calculate the Layer Output Matrix H
tempH = IW * P + BiasMatrix;

switch TF
    case 'sig'
        H = 1 ./ (1 + exp(-tempH));
    case 'sin'
        H = sin(tempH);
    case 'hardlim'
        H = hardlim(tempH);
    case 'softplus'
        H = log(1+exp(tempH));
end
% Calculate the Output Weight Matrix
LW = pinv(H') * T';
% c_rho = 2^27;
% out_w = pinv(speye(size(H,1))/c_rho + H*H') * (H*T');
% W_beta = 1 ./ max(abs(out_w),0.000001);
% W_beta = diag(sparse(W_beta));
% out_w = (W_beta/c_rho + H'*W*H) \ (H*(W*T'));
%% 
rankH=rank(H);
if rankH==size(H,1)
    %OutputWeight=inv(H * H'+diag(rand(size(H,1),1))) * H * T';              % output weights in case of full column rank
%     out_w=inv(H * H') * H * T';              % output weights in case of full column rank
    out_w=(H * H') \ H * T';              % output weights in case of full column rank
elseif rankH==size(H,2)
%     out_w= H *inv(H' * H) * T';              % output weights in case of full row rank
    out_w= H' *inv(H * H')\ T';              % output weights in case of full row rank
else 
        M = H * H';
        c = H * T';
        eps1 = 10^(-4);
        k = 1./size(M,1).* trace(M);
        B0 = k * eye(size(M));
        beta0 = 1./k * c;                           
        % or beta= inv(B0)*c
        if norm(c - B0 * beta0) > eps1
            b = B0 + (c - B0 * beta0) * (c - B0 * beta0)'./((c - B0 * beta0)'*beta0);       % rank-1 modification of B0
        else
            b = B0 + (c * c')./(c'* beta0) - (B0*(beta0*beta0')*B0)./(beta0'*B0*beta0);       % rank-2 modification of B0
        end
%         out_w = inv(B)*c;                   % output weights in case of not full rank
        out_w = b\c;                   % output weights in case of not full rank
end
%% 
    M0=inv(H*H');
    Beta0=M0*H*T';
    temp1M=M0*H*H'*M0;
    temp2M=H*H'*M0+1;
    M1=M0-temp1M/temp2M;
    tempBeta=M1*H*(T-(H'*Beta0)')';%%%%%%tempBeta «NumberofHiddenNeurons °¡j£¨
    Beta1=Beta0+tempBeta;
    Beta0=Beta1;