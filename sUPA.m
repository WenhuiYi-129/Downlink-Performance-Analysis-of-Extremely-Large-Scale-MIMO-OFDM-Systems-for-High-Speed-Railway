function [EDOF1,EDOF2]=sUPA(Nth,Ntv,Nrh,Nrv,L_tH,L_tV,L_rH,L_rV,fc,c,xr,yr,zr,u,l,s,f,L,Ts)

EDOF1=sUPAclosed(Nth,Ntv,Nrh,Nrv,L_tH,L_tV,L_rH,L_rV,fc,c,xr,yr,zr,u,l,s,f,L,Ts);
EDOF2=sUPAunclosed(Nth,Ntv,Nrh,Nrv,L_tH,L_tV,L_rH,L_rV,fc,c,xr,yr,zr,u,l,s,f,L,Ts);

end