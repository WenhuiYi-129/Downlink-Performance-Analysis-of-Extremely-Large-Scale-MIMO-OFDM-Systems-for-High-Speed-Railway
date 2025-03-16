function [H_ks]=channelmode_far(delta, UserSpacing, MH,MV,NH,NV,K_laisi,delta_phi,L,f_s,u,UsersNum,Ts,l,s,RecLength,TraLength)
%% 该函数生成基站到第k个用户的信道矩阵
%%输入参数设置依次为：UPA行数，UPA列数，水平方向角，俯仰角，路径损耗系数，第k个用户的距离，
%%参考距离，第n个簇的标准到达角，角度标准差，散射簇数量，Rician因子,发送载波索引，接收载波索引
%%s序列下的载波频率，矢量速度, 第k个用户的位置,符号持续时间
M=MH*MV;
N=NH*NV;
c=3e8;
lambda_s=c/f_s;
H_ks =zeros(N*UsersNum,M);

for user=1:UsersNum
    Lth=TraLength.L_x;
Ltv=TraLength.L_y;
Lrh=RecLength(user).L_x;
Lrv=RecLength(user).L_y;
    RecLength(user)=RecLength(1);
    RecLength(user).distance(1) =(user-1)*UserSpacing+RecLength(1).distance(1);
    distances=norm(RecLength(user).distance-TraLength.distance,2);
    [channelGaindB,~] = functionChannelgain(distances);
   channelGaindB=db2pow(channelGaindB);
    fd=Doppler(f_s,u,RecLength(user).distance,TraLength.distance,Ts);%%标准多普勒频移
    theta_k = asin((RecLength(user).distance(1)-TraLength.distance(1))/((RecLength(user).distance(3)-TraLength.distance(3))));
    phi_k = acos((RecLength(user).distance(2)-TraLength.distance(2))/distances);
    
    %% los径信道建模
    if (l-s+fd)==0
        I_k=1;
    else
    I_k=(sin(pi*(l-s+fd))/L/sin(pi*(l-s+fd)/L))*exp(1j*pi*(1-1/L)*(l-s+fd));
    end
    if l==s
    I_D=1;
    else
        I_D=(-1)^(l-s)*f_s*norm(u)/c/(sqrt(2)*(l-s));
    end
   
    h_sh=zeros(MH,1);
    h_sv=zeros(MV,1);
    h_rh=zeros(NH,1);
    h_rv=zeros(NV,1);
    for mh=1:MH
        h_sh(mh)=exp(-1j*pi*delta*sin(theta_k)*sin(phi_k)*(mh-1));
    end
    for mv=1:MV
        h_sv(mv)=exp(-1j*pi*delta*cos(phi_k)*(mv-1));
    end
    for rh=1:NH
        h_rh(rh)=exp(-1j*pi*delta*sin(theta_k)*sin(phi_k)*(rh-1));
    end
    for rv=1:NV
        h_rv(rv)=exp(-1j*pi*delta*cos(phi_k)*(rv-1));
    end
    vector_s=kron(h_sh,h_sv);
    vector_r=kron(h_rh,h_rv);
    H_los = sqrt(K_laisi/(K_laisi+1))*sqrt(channelGaindB)*vector_r*vector_s';
  
  %% nlos径信道建模
  path_num = 6;
 h_nlos = zeros(N,M);

theta_s = rand(path_num,1)*pi;                                          % \theta_t
phi_s   = rand(path_num,1)*pi;                                          % \phi_t
theta_r = rand(path_num,1)*pi;                                          % \theta_r
phi_r   = rand(path_num,1)*pi;                                          % \phi_r
z = (randn(1, path_num) + 1i * randn(1,path_num));  % 随机生成复数数列，实部和虚部均为正态分布
% 计算当前数列的2范数
current_norm = norm(z, 'fro');

path_loss= sqrt(1/(1+K_laisi))*sqrt(channelGaindB)*z*(1/current_norm); % 生成复数序列，模为1

A_s = zeros(M,1);
A_r = zeros(1,N);
for th=1:path_num
k_s = (2*pi/lambda_s)*[cos(theta_s(th))*cos(phi_s(th)), cos(theta_s(th))*sin(phi_s(th)), sin(theta_s(th))];
k_r = (2*pi/lambda_s)*[cos(theta_r(th))*cos(phi_r(th)), cos(theta_r(th))*sin(phi_r(th)), sin(theta_r(th))];
for m=1:M
    i=mod(m-1,MH);
    j=floor((m-1)/MH);
    r_tm=TraLength.distance+[-Lth/2+i*delta -Ltv/2+j*delta 0];
    A_s(m,1)=exp(1j*dot(k_s,r_tm));
end

for n=1:N
    t=mod(n-1,NH);
    k=floor((n-1)/NH);
    r_rn=RecLength(user).distance+[-Lrh/2+t*delta -Lrv/2+k*delta 0];
    A_r(1,n) = exp(-1j*dot(k_r,r_rn));
end

  h_nlos= h_nlos+sqrt(1/path_num)*path_loss(th)*A_r*A_s;

% path_loss= sqrt(1/(1+K_laisi))*sqrt(channelGaindB); % 生成复数序列，模为1
% h_nlos_1=zeros(N,M);
% for th = 1:path_num
%  h_sh_2=zeros(MH,1);
%     h_sv_2=zeros(MV,1);
%     h_rh_2=zeros(NH,1);
%     h_rv_2=zeros(NV,1);
%     for mh=1:MH
%         h_sh_2(mh)=exp(-1j*pi*delta*sin(theta_t(th))*sin(phi_t(th))*(mh-1));
%     end
%     for mv=1:MV
%         h_sv_2(mv)=exp(-1j*pi*delta*cos(phi_t(th))*(mv-1));
%     end
%     for rh=1:NH
%         h_rh_2(rh)=exp(-1j*pi*delta*sin(theta_r(th))*sin(phi_r(th))*(rh-1));
%     end
%     for rv=1:NV
%         h_rv_2(rv)=exp(-1j*pi*delta*cos(phi_r(th))*(rv-1));
%     end
%      vector_s_2=kron(h_sh_2,h_sv_2);
%     vector_r_2=kron(h_rh_2,h_rv_2);
%     h_nlos_1=h_nlos_1+(1/sqrt(path_num))*vector_r_2*vector_s_2';
% end
% h_nlos = h_nlos + path_loss*h_nlos_1/100;
end
    H_ks((user-1)*N+1:user*N,:)= I_k*H_los +I_D*h_nlos;
end

end
