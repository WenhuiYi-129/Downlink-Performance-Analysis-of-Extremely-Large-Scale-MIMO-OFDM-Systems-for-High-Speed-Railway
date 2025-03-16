function [H_channel]=channel(f,transmit_position,receive_position,Lth,Ltv,Lrh,Lrv,Mh,Mv,Nh,Nv,l,s,u,L,Ts)
c=3e8;
M=Mh*Mv;
N=Nh*Nv;
H_channel=zeros(N,M);
deltath=Lth/Mh;
deltatv=Ltv/Mv;
deltarh=Lrh/Nh;
deltarv=Lrv/Nv;
for n=1:N
    t=mod(n-1,Nh);
    k=floor((n-1)/Nh);
    for m=1:M
        i=mod(m-1,Mh);
        j=floor((m-1)/Mh);
        rtm=transmit_position+[0 -Lth/2+i*deltath -Ltv/2+j*deltatv];
        rrn=receive_position+[0 -Lrh/2+t*deltarh -Lrv/2+k*deltarv];
        fd=Doppler(f,u,rrn,rtm,Ts);
        d=norm(rrn-rtm);
        h=exp(-1j*2*pi*f*d/c)/(4*pi)/d;
        if l-s+fd~=0
        h=h*(sin(pi*(l-s+fd))/L/sin(pi*(l-s+fd)/L))*exp(1j*pi*(1-1/L)*(l-s+fd));
        end
            
        H_channel(n,m)=h;        
    end
end

end