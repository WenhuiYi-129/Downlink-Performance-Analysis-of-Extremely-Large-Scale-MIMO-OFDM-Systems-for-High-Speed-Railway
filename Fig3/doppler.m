function [w,ending]=doppler(rr,rt,u,fc)
c=3e8;
if norm(u)~=0
phi=dot((rr-rt),u)/norm(rr-rt)/norm(u);
ending=-fc*norm(u)*phi/c;
else
    ending=0;
end
w=fc*norm(u)/c;