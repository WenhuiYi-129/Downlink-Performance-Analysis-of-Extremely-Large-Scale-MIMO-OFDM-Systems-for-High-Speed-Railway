clc
clear all
tic
global FRAME
FRAME=20;
double=0; % the double scattering environment 
step=2;
% SNR_Range=-10:step:20;

% global spacing % the spacing between RIS elements
global lambda
global RecLength
global TraLength
fc=2.1e9;
c=3e8;
lambda=c/fc; % wavelength 
Ts=67e-6;
L=8;
power=0.2;
for l=1:L
    f(l)=fc+(2*l-L-1)/2/Ts;
end
M=[8 12 16 20];
lenM=length(M);

D=4000;
dtransmit_position=[0 0 0];
itransmit1_position=[-D 0 0];
itransmit2_position=[D 0 0];
v1=0:100:600;
lenv=length(v1);
UsersNum=1; 
UserSpacing=100*lambda;
SINRend1=zeros(lenv,length(M));
SEend1=zeros(lenv,length(M));
for i=1:lenM
    disp(M(i));
Nr_X=M(i); Nr_Y=M(i);
Ns_X=M(i); Ns_Y=M(i);
RecSpacing=lambda/2;
TraSpacing=lambda/2;
RecLength.L_x=Nr_X*RecSpacing;RecLength.L_y=Nr_Y*RecSpacing;
RecLength.distance=[-1000,-6,2];
TraLength.L_x=Ns_X*TraSpacing;TraLength.L_y=Ns_Y*TraSpacing;
TraLength.distance=[0,0,0];
d=norm(RecLength.distance);
for j=1:lenv
    v=v1(j);
    disp(v);
       v=v*1000/3600;
    u=[v 0 0];
% laisi_k=13-0.03*d;
laisi_k=db2pow(30);
% laisi_k=100;
%% Generate the channel matrix
% Obtain the sampled eclipse and their coordinates given a surface
[RecSamplePoint,WD_RecX_vec,WD_RecY_vec]=WD_eclipse_sample(lambda,RecLength);
[TraSamplePoint,WD_TraX_vec,WD_TraY_vec]=WD_eclipse_sample(lambda,TraLength);

[RecVarianceVec,WD_RecX_vec,WD_RecY_vec,nr_act]=WaveDomainChannel(RecSamplePoint,WD_RecX_vec,WD_RecY_vec,lambda,RecLength);
[TraVarianceVec,WD_TraX_vec,WD_TraY_vec,ns_act]=WaveDomainChannel(TraSamplePoint,WD_TraX_vec,WD_TraY_vec,lambda,TraLength);

% average
if double~=1
    tt1=mean(TraVarianceVec);
    TraVarianceVec=tt1*ones(ns_act,1);
end
ChannelSigma=RecVarianceVec*TraVarianceVec';
[TraResVector,TraNumNs]=ResponseGenerate(TraSpacing,TraSpacing,lambda,TraLength,WD_TraX_vec,WD_TraY_vec,'Transmit'); 


% Multi-user Case
RecResVector=[];
for user=1:UsersNum 
    RecLength(user)=RecLength(1);
    RecLength(user).distance(1) =(user-1)*UserSpacing+RecLength(1).distance(1);
    tempUserInfo=RecLength(user);
    [RecResVectorPart,RecNumNr]=ResponseGenerate(RecSpacing,RecSpacing,lambda,tempUserInfo,WD_RecX_vec,WD_RecY_vec,'Receive');
    RecResVector=blkdiag(RecResVector,RecResVectorPart);
end

%% Pass Channel 
for frame=1:FRAME
TransSig=randn(UsersNum*RecNumNr,1)+1i*randn(UsersNum*RecNumNr,1); % transmitted signal

        distances=d;
        [channelGaindB,~] = functionChannelgain(distances);
        channelGain=db2pow(channelGaindB);
        ChannelRandom=sqrt(1/2)*(randn(UsersNum*nr_act,ns_act)+1i*randn(UsersNum*nr_act,ns_act));
              
        ChannelHa=M(i)^2*sqrt(channelGain)*ChannelSigma.*ChannelRandom/sum(ChannelSigma(:));
        Channel= sqrt(1/(1+laisi_k))*RecResVector*ChannelHa*TraResVector;
    [SINRend,SEend]=SINR(fc,Ts,dtransmit_position,itransmit1_position,itransmit2_position,RecLength.distance,TraLength.L_x,TraLength.L_y,RecLength.L_x,RecLength.L_y,M(i),M(i),M(i),M(i),u,L,power,Channel,laisi_k);
SINRend1(j,i)=SINRend1(j,i)+SINRend;
SEend1(j,i)=SEend1(j,i)+SEend;
end
SINRend1(j,i)=SINRend1(j,i)/FRAME;
SEend1(j,i)=SEend1(j,i)/FRAME;
 v2(j)=Doppler(fc,u,[-2000 -6 2],[0 0 0],Ts);
end 
end
figure(1)
h1=plot(v1,SINRend1(:,1),'-','LineWidth',2);
 hold on;  
  h2=plot(v1,SINRend1(:,2),'-','LineWidth',2);
 hold on; 
 h3=plot(v1,SINRend1(:,3),'-','LineWidth',2);
 hold on;  
  h4=plot(v1,SINRend1(:,4),'-','LineWidth',2);
 hold on; 
xlabel('normalized DFO','Interpreter','Latex');
ylabel('SINR $(\mathrm{dB})$','Interpreter','Latex');
legend([h1 h2 h3 h4],{'$M_H\times M_V=8\times 8$','$M_H\times M_V=12\times 12$','$M_H\times M_V=16\times 16$','$M_H\times M_V=20\times 20$'},'Interpreter','Latex','Location','Northwest')
set(gca,'FontSize',12);
grid on;
toc;
disp(['运行时间:',num2str(toc)]);
% figure(2)
% h1=plot(v1,SEend1(:,1),'-','LineWidth',2);
%  hold on;  
%   h2=plot(v1,SEend1(:,2),'-','LineWidth',2);
%  hold on; 
% xlabel('normalized DFO','Interpreter','Latex');
% ylabel('SINR $(\mathrm{dB})$','Interpreter','Latex');
% legend([h1 h2 h3 h4],{'$M_H\times M_V=2\times 2$','$M_H\times M_V=8\times 8$','$M_H\times M_V=14\times 14$','$M_H\times M_V=20\times 20$'},'Interpreter','Latex','Location','Northwest')
% set(gca,'FontSize',12);
% grid on;