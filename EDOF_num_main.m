clc
clear all
tic
global FRAME
FRAME=200;
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
for l=1:L
    f(l)=fc+(2*l-L-1)/2/Ts;
end
M=1:1:5;
lenM=length(M);
LT1=[25*lambda 31*lambda];
% LT1=[2*lambda 45*lambda];
v=350;
v=v*1000/3600;
u=[v 0 0];
UsersNum=1; 
UserSpacing=100*lambda;
EDOF1_sUPAclosed=zeros(length(LT1),lenM);
EDOF1_sUPAunclosed=zeros(length(LT1),lenM);
EDOF1_sULAclosed=zeros(length(LT1),lenM);
EDOF1_sULAunclosed=zeros(length(LT1),lenM);
for j=1:length(LT1)

for i=1:lenM
    disp(M(i));
Nr_X=M(i); Nr_Y=M(i);
Ns_X=M(i); Ns_Y=M(i);
    RecSpacing=LT1(j)*sqrt(2)/2/M(i);
    TraSpacing=LT1(j)*sqrt(2)/2/M(i);
RecLength.L_x=LT1(j)*sqrt(2)/2;RecLength.L_y=LT1(j)*sqrt(2)/2;
RecLength.distance=[20*lambda,-6,2];

TraLength.L_x=LT1(j)*sqrt(2)/2;TraLength.L_y=LT1(j)*sqrt(2)/2;
TraLength.distance=[0,0,0];
d=norm(RecLength.distance);

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

TransSig=randn(UsersNum*RecNumNr,1)+1i*randn(UsersNum*RecNumNr,1); % transmitted signal

        distances=d;
        [channelGaindB,~] = functionChannelgain(distances);
        channelGain=db2pow(channelGaindB);
%     for frame=1:FRAME
%         noise=sqrt(var_noise/2)*(randn(UsersNum*RecNumNr,1)+1i*randn(UsersNum*RecNumNr,1));
        ChannelRandom=sqrt(1/2)*(randn(UsersNum*nr_act,ns_act)+1i*randn(UsersNum*nr_act,ns_act));
              
        ChannelHa=M(i)^2*sqrt(channelGain)*ChannelSigma.*ChannelRandom/sum(ChannelSigma(:));
%         ChannelHa= diag(RecVarianceVec)*ChannelRandom*diag(TraVarianceVec);
%         Channel=1/sqrt(UsersNum*RecNumNr*TraNumNs)*RecResVector*ChannelHa*TraResVector;
        Channel= sqrt(1/(1+laisi_k))*RecResVector*ChannelHa*TraResVector;
%                 Channel= RecResVector*ChannelHa*TraResVector;

l=2;
s=2;
EDOF1_sUPAclosed(j,i)=sUPAclosed(Ns_X,Ns_Y,Nr_X,Nr_Y,TraLength.L_x,TraLength.L_y,RecLength.L_x,RecLength.L_y,fc,c,RecLength.distance(1),RecLength.distance(2),RecLength.distance(3),u,l,s,f,L,Ts,Channel,laisi_k,ChannelSigma,TraResVector');
EDOF1_sUPAunclosed(j,i)=sUPAunclosed(Ns_X,Ns_Y,Nr_X,Nr_Y,TraLength.L_x,TraLength.L_y,RecLength.L_x,RecLength.L_y,fc,c,RecLength.distance(1),RecLength.distance(2),RecLength.distance(3),u,l,s,f,L,Ts,Channel,laisi_k);
        
 
end 
end

% %%ULA
% for j=1:length(LT1)
% 
% for i=1:lenM
% Nr_X=1; Nr_Y=M(i)^2;
% Ns_X=1; Ns_Y=M(i)^2;
% spacing_x=lambda/2;
%     RecSpacing=LT1(j)/M(i)^2;
%     TraSpacing=LT1(j)/M(i)^2;
% RecLength.L_x=Nr_X*spacing_x;RecLength.L_y=Nr_Y*RecSpacing;
% RecLength.distance=[20*lambda,-6,2];
% 
% TraLength.L_x=Ns_X*spacing_x;TraLength.L_y=Ns_Y*TraSpacing;
% TraLength.distance=[0,0,0];
% d=norm(RecLength.distance);
% % laisi_k=13-0.03*d;
% % laisi_k=db2pow(20);
% laisi_k=db2pow(30);
% %% Generate the channel matrix
% % Obtain the sampled eclipse and their coordinates given a surface
% [RecSamplePoint,WD_RecX_vec,WD_RecY_vec]=WD_eclipse_sample(lambda,RecLength);
% [TraSamplePoint,WD_TraX_vec,WD_TraY_vec]=WD_eclipse_sample(lambda,TraLength);
% 
% [RecVarianceVec,WD_RecX_vec,WD_RecY_vec,nr_act]=WaveDomainChannel(RecSamplePoint,WD_RecX_vec,WD_RecY_vec,lambda,RecLength);
% [TraVarianceVec,WD_TraX_vec,WD_TraY_vec,ns_act]=WaveDomainChannel(TraSamplePoint,WD_TraX_vec,WD_TraY_vec,lambda,TraLength);
% 
% % average
% if double~=1
%     tt1=mean(TraVarianceVec);
%     TraVarianceVec=tt1*ones(ns_act,1);
% end
% ChannelSigma=RecVarianceVec*TraVarianceVec';
% [TraResVector,TraNumNs]=ResponseGenerate(spacing_x,TraSpacing,lambda,TraLength,WD_TraX_vec,WD_TraY_vec,'Transmit'); 
% 
% 
% % Multi-user Case
% RecResVector=[];
% for user=1:UsersNum 
%     RecLength(user)=RecLength(1);
%     RecLength(user).distance(1) =(user-1)*UserSpacing+RecLength(1).distance(1);
%     tempUserInfo=RecLength(user);
%     [RecResVectorPart,RecNumNr]=ResponseGenerate(spacing_x,RecSpacing,lambda,tempUserInfo,WD_RecX_vec,WD_RecY_vec,'Receive');
%     RecResVector=blkdiag(RecResVector,RecResVectorPart);
% end
% 
% %% Pass Channel 
% 
% TransSig=randn(UsersNum*RecNumNr,1)+1i*randn(UsersNum*RecNumNr,1); % transmitted signal
% 
%         distances=d;
%         [channelGaindB,~] = functionChannelgain(distances);
%         channelGain=db2pow(channelGaindB);
% %     for frame=1:FRAME
% %         noise=sqrt(var_noise/2)*(randn(UsersNum*RecNumNr,1)+1i*randn(UsersNum*RecNumNr,1));
%         ChannelRandom=sqrt(1/2)*(randn(UsersNum*nr_act,ns_act)+1i*randn(UsersNum*nr_act,ns_act));
%               
% %         ChannelHa=ChannelSigma.*ChannelRandom;
% %         ChannelHa= diag(RecVarianceVec)*ChannelRandom*diag(TraVarianceVec);
%                ChannelHa=M(i)^2*sqrt(channelGain)*ChannelSigma.*ChannelRandom/sum(ChannelSigma(:));
% %         Channel=1/sqrt(UsersNum*RecNumNr*TraNumNs)*RecResVector*ChannelHa*TraResVector;
%         Channel= sqrt(1/(1+laisi_k))*RecResVector*ChannelHa*TraResVector;
%  l=2;
% s=2;
% EDOF1_sULAclosed(j,i)=sULAclosed(Ns_Y,Nr_Y,TraLength.L_y,RecLength.L_y,fc,c,RecLength.distance(1),RecLength.distance(2),RecLength.distance(3),u,l,s,f,L,Ts,Channel,laisi_k,ChannelSigma);
% EDOF1_sULAunclosed(j,i)=sULAunclosed(Ns_Y,Nr_Y,TraLength.L_y,RecLength.L_y,fc,c,RecLength.distance(1),RecLength.distance(2),RecLength.distance(3),u,l,s,f,L,Ts,Channel,laisi_k);
%      
% end 
% end

figure(1);
h1=plot(M.^2,EDOF1_sUPAclosed(1,:),'o','LineWidth',1.5);
hold on ;
h2=plot(M.^2,EDOF1_sUPAunclosed(1,:),'--','LineWidth',2);
hold on ;

h3=plot(M.^2,EDOF1_sUPAclosed(2,:),'o','LineWidth',1.5);
hold on ;
h4=plot(M.^2,EDOF1_sUPAunclosed(2,:),'--','LineWidth',2);
hold on ;
h5=plot(M.^2,EDOF1_sULAclosed(1,:),'o','LineWidth',1.5);
hold on ;
h6=plot(M.^2,EDOF1_sULAunclosed(1,:),'--','LineWidth',2);
hold on ;

h7=plot(M.^2,EDOF1_sULAclosed(2,:),'o','LineWidth',1.5);
hold on ;
h8=plot(M.^2,EDOF1_sULAunclosed(2,:),'--','LineWidth',2);
hold on ;
xlabel('antennes number','Interpreter','Latex');
ylabel('EDoF','Interpreter','Latex');
legend([h6 h2 h8 h4 h1],{'ULA,$D=25\lambda$','UPA,$D=25\lambda$','ULA,$D=35\lambda$','UPA,$D=35\lambda$','Analytical'},'Interpreter','Latex','Location','Northwest')
set(gca,'FontSize',12);
grid on;
toc;
disp(['运行时间:',num2str(toc)]);
  
