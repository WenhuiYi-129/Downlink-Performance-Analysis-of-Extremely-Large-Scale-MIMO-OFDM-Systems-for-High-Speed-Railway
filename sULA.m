function [EDOF1,EDOF2]=sULA(Nt,Nr,Lt,Lr,fc,c,xr,yr,zr,u,l,s,f,L,Ts)

EDOF1=sULAclosed(Nt,Nr,Lt,Lr,fc,c,xr,yr,zr,u,l,s,f,L,Ts);
EDOF2=sULAunclosed(Nt,Nr,Lt,Lr,fc,c,xr,yr,zr,u,l,s,f,L,Ts);

end