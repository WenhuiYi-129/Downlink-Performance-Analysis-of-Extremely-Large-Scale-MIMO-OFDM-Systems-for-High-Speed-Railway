function [G]=greenfunction(rr,rt,f,flag)
%%Flag represents different polarization modes (Scalar,single polarization, dual polarization, tertiary)
c=3e8;
wavenumber=2*pi*f/c;
d=norm(rr-rt);
G=exp(-1j*wavenumber*d)/(4*pi*d);
  ar=(rr-rt)./d;   %Unit direction vector
  I=eye(3,3);
  G1=1+1j/(wavenumber*d)-1/(wavenumber*d)^2;
  G2=3/(wavenumber*d)^2-3j/(wavenumber*d)-1;
 Gd=G1*I*G+G2*(ar'*ar)*G;
 G=Gd(1:flag,1:flag);

end