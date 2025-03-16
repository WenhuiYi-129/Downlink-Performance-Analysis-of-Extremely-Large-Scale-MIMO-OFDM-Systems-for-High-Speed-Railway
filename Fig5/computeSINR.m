function [SINRs2]=computeSINR(fc,Ts,dtransmit_position,itransmit1_position,itransmit2_position,receive_position,Lth,Ltv,Lrh,Lrv,Mh,Mv,Nh,Nv,s,u,L,power,Channel,laisi_k)
% P_=sqrt(power)*eye(Mh*Mv);
%  p_power=norm(P_,'fro')^2;
 %Noise figure
  noisefigure=9;
  %communication bandwidth
  B=20e6;
  noisepowerdbm=-174+10*log10(B)+noisefigure;%%noise power
  noisepower=db2pow(noisepowerdbm)/1000;
%   power=0.2;
% noisepower=norm(noisepower,'fro')^2*L;
N=Nh*Nv;
  dpower=0;
  ipower=0;
  noisepower=noisepower*N*L;
for l=1:L
    f=fc+(2*l-L-1)/2/Ts;
   H_channel=channel(f,dtransmit_position,receive_position,Lth,Ltv,Lrh,Lrv,Mh,Mv,Nh,Nv,l,s,u,L,Ts,Channel,laisi_k);
    if l==s
dpower=dpower+norm(sqrt(power)*H_channel,'fro')^2;
    else
         ipower=ipower+norm(sqrt(power)*H_channel,'fro')^2;
    end
end
SINRs2=dpower/(ipower+noisepower);
end
