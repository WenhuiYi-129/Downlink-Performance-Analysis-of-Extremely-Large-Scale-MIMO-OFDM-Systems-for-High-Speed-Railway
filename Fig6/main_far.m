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
delta_phi=5*pi/180;
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

 [RecLength,ChannelSigma,RecResVector,TraResVector,nr_act,ns_act]=channel_generate(laisi_k,UserSpacing,RecSpacing,TraSpacing,lambda,RecLength,TraLength,UsersNum,double);

 

power=0.2;


for v=v1(1):100:v1(end)
    disp(v);
    var_noise=noisepower;
    v2=v*1000/3600;
    v2=0;
    u=[v2 0 0];
    
       MRT_collect=[];
    ZF_collect=[];
    MMSE_collect=[];
    ZF_NS_collect=[];
    
    MRT_theo_collect=[];
    ZF_iter_collect=[];
    MMSE_theo_collect=[];
    alpha_mrt_collect=[];
for frame=1:FRAME
M=Ns_X*Ns_Y;
 N=Nr_X*Nr_Y;  
H_channel=zeros(N*UsersNum,M,L,L);
b=toc;
fprintf('frame=%d, speed=%f, times = %d\n', frame, v, b);
for l =1:L
    for s=1:L
 
 H_TOTAOL=channelmode_far(RecSpacing,UserSpacing,Ns_X,Ns_Y,Nr_X,Nr_Y,laisi_k,delta_phi,L,f1(s),u,UsersNum,Ts,l,s,RecLength,TraLength);
    H_channel(:,:,l,s)= H_TOTAOL;
    end
end
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
figure(1)
p1=plot(v1,abs(MRT_snr),'s-m','linewidth',2);
hold on
p2=plot(v1,abs(ZF_snr),'o-b','linewidth',2);
hold on 
p3=plot(v1,abs(MMSE_snr),'d-r','linewidth',2);
hold on
legend([p1 p2 p3],'MRT(16 $\times$ 16) ','ZF(16 $\times$ 16) ','MMSE(16 $\times$ 16)','MRT(12 $\times$ 12) ','ZF(12 $\times$ 12) ','MMSE(12 $\times$ 12)','Interpreter','latex' )

xlabel('speed $v$ (km/h)','Interpreter','latex')
ylabel('Spectral Efficiency (bits/s/Hz)','Interpreter','latex')
grid on

  
