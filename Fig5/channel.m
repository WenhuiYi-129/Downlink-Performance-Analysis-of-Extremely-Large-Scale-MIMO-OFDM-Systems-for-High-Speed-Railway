function [H_channel]=channel(f,transmit_position,receive_position,Lth,Ltv,Lrh,Lrv,Mh,Mv,Nh,Nv,l,s,u,L,Ts,Channel,laisi_k)
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
        rtm=transmit_position+[-Lth/2+i*deltath -Ltv/2+j*deltatv 0];
        rrn=receive_position+[-Lrh/2+t*deltarh -Lrv/2+k*deltarv 0];
        fd=Doppler(f,u,rrn,rtm,Ts);
        d=norm(rrn-rtm);
     [channelGaindB,~] = functionChannelgain(d);
        channelGain=sqrt(db2pow(channelGaindB));
        h=exp(-1j*2*pi*f*d/c)*channelGain*sqrt(laisi_k/(1+laisi_k));
        if l-s+fd~=0
        h=h*(sin(pi*(l-s+fd))/L/sin(pi*(l-s+fd)/L))*exp(1j*pi*(1-1/L)*(l-s+fd));
        end
        w=f*norm(u)*Ts/c;
         if l~=s
               Channel(n,m)=Channel(n,m)*(-1)^(l-s)*w/sqrt(2)/(l-s);
          end   
        H_channel(n,m)=h;        
    end
end
H_channel=H_channel+Channel;
end