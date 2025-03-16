function [SINR1,SEs1]=SINR(fc,Ts,dtransmit_position,itransmit1_position,itransmit2_position,receive_position,Lth,Ltv,Lrh,Lrv,Mh,Mv,Nh,Nv,u,L,power,Channel,laisi_k)
SINR1=0;
SINR2=0;
SEs1=0;
SEs2=0;
for s=1:L
    [sinrm1]=computeSINR(fc,Ts,dtransmit_position,itransmit1_position,itransmit2_position,receive_position,Lth,Ltv,Lrh,Lrv,Mh,Mv,Nh,Nv,s,u,L,power,Channel,laisi_k);
    SEs1=SEs1+log2(1+sinrm1);
    SINR1=SINR1+sinrm1;
end
SINR1=SINR1/L;
SINR1=pow2db(SINR1);
end