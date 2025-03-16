clc
clear all
tic
global FRAME
FRAME=300;
double=0; % the double scattering environment 
step=2;
%Noise figure
  noisefigure=9;
  %communication bandwidth
  B=20e6;
  noisepowerdbm=-174+10*log10(B)+noisefigure;%%noise power
  noisepower=db2pow(noisepowerdbm)/1000;
% SNR_Range=18:step:20;
% power1=10:20:100;
% power1=power1./1000;
v1=0:100:600;
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
    f1(l)=fc+(2*l-L-1)/2/Ts;
end
UsersNum=4; 
UserSpacing=25;
RecSpacing_num=2;
RecSpacing=lambda/RecSpacing_num;
TraSpacing_num=2;
TraSpacing=lambda/TraSpacing_num;
laisi_k=db2pow(30);
Nr_X=10; Nr_Y=10;
Ns_X=10; Ns_Y=10;
% v=350*1000/3600;
% u=[v 0 0];
RecLength.L_x=Nr_X*RecSpacing;RecLength.L_y=Nr_Y*RecSpacing;
RecLength.distance=[-100,-6,2];

TraLength.L_x=Ns_X*TraSpacing;TraLength.L_y=Ns_Y*TraSpacing;
TraLength.distance=[0,0,0];
% % The length of received  surface
% RecLength =struct('L_x',2*lambda,'L_y',4*lambda,'distance',[6,4,2]);
% % The length of transmitted  surface
% TraLength=struct('L_x',16*lambda,'L_y',10*lambda,'distance',[0,0,0]);
%% Generate the channel matrix
[RecLength,ChannelSigma,RecResVector,TraResVector,nr_act,ns_act]=channel_generate(laisi_k,UserSpacing,RecSpacing,TraSpacing,lambda,RecLength,TraLength,UsersNum,double);

MRT_snr=[];
ZF_snr=[];
ZF_NS_snr_theo=[];
MMSE_snr=[];
MRT_snr_theo=[];
ZF_snr_theo=[];
MMSE_snr_theo=[];
% for SNR=SNR_Range(1):step:SNR_Range(end)
% for power=power1(1):0.02:power1(end)
% power=0.1;
% P=50;

%     SNR1=db2pow(SNR);
%     power=SNR1*var_noise;

 Channel=zeros(Nr_X*Nr_Y*UsersNum,Ns_X*Ns_Y,FRAME);
    for frame=1:FRAME
          ChannelRandom=sqrt(1/2)*(randn(UsersNum*nr_act,ns_act)+1i*randn(UsersNum*nr_act,ns_act));
      ChannelHa=Ns_X^2*ChannelSigma.*ChannelRandom/sum(ChannelSigma(:));
        Channel(:,:,frame)=  RecResVector*ChannelHa*TraResVector;
    end
disp('NLOS end');
toc;
disp(['运行时间:',num2str(toc)]);    
 
    
power=0.2;
for v=v1(1):100:v1(end)
    v=0;
    disp(v);
    var_noise=noisepower;
    v=v*1000/3600;
    u=[v 0 0];
    
       MRT_collect=[];
    ZF_collect=[];
    MMSE_collect=[];
    ZF_NS_collect=[];
    
    MRT_theo_collect=[];
    ZF_iter_collect=[];
    MMSE_theo_collect=[];
    alpha_mrt_collect=[];
for frame=1:FRAME
N=Nr_X*Nr_Y;
M=Ns_X*Ns_Y;
H_channel=zeros(Nr_X*Nr_Y*UsersNum,Ns_X*Ns_Y,L,L);
for s=1:L
    transmit_position=TraLength.distance;
    for l=1:L
        LOSH=[];
        H=zeros(Nr_X*Nr_Y,Ns_X*Ns_Y);
    for user=1:UsersNum 
    receive_position=RecLength(user).distance;      
     H=channel(f1(s),transmit_position,receive_position,TraLength.L_x,TraLength.L_y,RecLength(user).L_x,RecLength(user).L_y,Ns_X,Ns_Y,Nr_X,Nr_Y,l,s,u,L,Ts,0,laisi_k);
      LOSH=[LOSH,transpose(H)];
    end
   
    LOSH=transpose(LOSH);
    w=f1(s)*norm(u)*Ts/c;
         if l~=s
               Channel_nlos=Channel(:,:,frame)*(-1)^(l-s)*w/sqrt(2)/(l-s);
         else
             Channel_nlos=Channel(:,:,frame);
         end   
       LOSH=LOSH+Channel_nlos;
       H_channel(:,:,s,l)=LOSH;
   end
end 
    disp(['speed=',num2str(v)]);
    disp(['frame=',num2str(frame)]);
toc;
disp(['运行时间:',num2str(toc)]); 
        % System Model
        [BFMatrix_MMSE,MMSE_iter]=BFSchemes_MMSE (H_channel,var_noise,power,UsersNum,L,Nr_X*Nr_Y,Ns_X*Ns_Y);
        MMSE_collect=[MMSE_collect,MMSE_iter];       
        [BFMatrix_ZF,ZF_iter]=BFSchemes_ZF (H_channel,var_noise,power,UsersNum,L,Nr_X*Nr_Y,Ns_X*Ns_Y);
        ZF_collect=[ZF_collect,ZF_iter];      
        [BFMatrix_MRT,MRT_iter]=BFSchemes_MRT(H_channel,var_noise,power,UsersNum,L,Nr_X*Nr_Y,Ns_X*Ns_Y);
        MRT_collect=[MRT_collect,MRT_iter];      
end
    MRT_snr=[MRT_snr,sum(MRT_collect)/FRAME];
    ZF_snr=[ZF_snr,sum(ZF_collect)/FRAME];
    MMSE_snr=[MMSE_snr,sum(MMSE_collect)/FRAME];  
end
toc;
disp(['运行时间:',num2str(toc)]);
% 
figure
p1=plot(v1,abs(MRT_snr),'s-m','linewidth',2);
hold on
p2=plot(v1,abs(ZF_snr),'o-b','linewidth',2);
hold on 
p3=plot(v1,abs(MMSE_snr),'d-r','linewidth',2);
hold on
legend([p1 p2 p3],'MRT','ZF ','MMSE','Interpreter','latex' )

xlabel('speed $v$ (km/h)','Interpreter','latex')
ylabel('Spectral Efficiency (bits/s/Hz)','Interpreter','latex')
grid on

  
