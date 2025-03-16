function [EDOF]=sULAunclosed(Nt,Nr,Lt,Lr,fc,c,xr,yr,zr,u,l,s,f,L,Ts,Channel,laisi_k)

%%unclosed  
     M=Nt;
     N=Nr;
        intervalt=Lt/M;
        intervalr=Lr/N;
        H=zeros(N,M);
                distances=sqrt(xr^2+yr^2+zr^2);
        [channelGaindB,~] = functionChannelgain(distances);
        channelGain=sqrt(db2pow(channelGaindB));
%         fd=doppler(location(0,0,0),location(xr,yr,zr),u,fc)*Ts;
        for n=1:N
            for m =1:M
               R_t_m=location(0,-Lt/2+(m-1)*intervalt,0);
               R_r_n=location(xr,yr-Lr/2+(n-1)*intervalr,zr);
                [w,fd]=doppler(R_t_m,R_r_n,u,f(l));
               w=w*Ts;
               fd=fd*Ts;
               if l~=s
                   Channel(n,m)=Channel(n,m)*(-1)^(l-s)*w/sqrt(2)/(l-s);
               end     
               k0=2*pi*f(l)/c;
               D=norm(R_r_n-R_t_m);
               if (l-s+fd)==0
        G=sqrt(laisi_k/(1+laisi_k))*complex(cos(-k0*D),sin(-k0*D))*channelGain*(distances/D);
               else
                G=sqrt(laisi_k/(1+laisi_k))*complex(cos(-k0*D),sin(-k0*D))*channelGain*(distances/D)*DFO(l,s,fd,L);
               end
               H(n,m)=G;
            end
        end
         H=H+Channel;
        R=H'*H;
        numerator=0;
        for ff=1:M
            numerator=numerator+R(ff,ff);
        end
        denominator=0;
        for ttt=1:M
            for ff=1:M
                denominator=denominator+abs(R(ttt,ff))^2;
            end
        end
        EDOF=(abs(numerator)^2)/denominator;

end