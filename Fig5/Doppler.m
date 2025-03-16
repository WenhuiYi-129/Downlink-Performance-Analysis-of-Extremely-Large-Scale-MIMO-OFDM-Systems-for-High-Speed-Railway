function [ending]=Doppler(f,u,rr,rt,Ts)
v=norm(u);
c=3e8;
if v==0
    ending=0;
else
    phi=dot(rr-rt,u)/norm(rr-rt)/v;
    ending=-f*v*phi*Ts/c;
    
end
end