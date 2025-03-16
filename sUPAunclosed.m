function [EDOF]=sUPAunclosed(Nth,Ntv,Nrh,Nrv,L_tH,L_tV,L_rH,L_rV,fc,c,xr,yr,zr,u,l,s,f,L,Ts,Channel,laisi_k)
%%unclosed  
      TX_number_h=Nth;
        TX_number_r=Ntv;
        RX_number_h=Nrh;
        RX_number_r=Nrv;
        M=TX_number_h*TX_number_r;%%Number of transmitter antennas
        N=RX_number_h*RX_number_r;%%Number of receiver antennas
        delta_th=L_tH/TX_number_h;
        delta_tv=L_tV/TX_number_r;
        delta_rh=L_rH/RX_number_h;
        delta_rv=L_rV/RX_number_r;
        H=zeros(N,M);
        distances=sqrt(xr^2+yr^2+zr^2);
        [channelGaindB,~] = functionChannelgain(distances);
        channelGain=sqrt(db2pow(channelGaindB));
        for m=1:M
            i1=mod(m-1,TX_number_h);
                j1=floor((m-1)/TX_number_h);
            for n =1:N
                 t1=mod(n-1,RX_number_h);
                 k1=floor((n-1)/RX_number_h);
               R_t_m=location(-L_tH/2+i1*delta_th,-L_tV/2+j1*delta_tv,0);
               R_r_n=location(xr-L_rH/2+t1*delta_rh,yr-L_rV/2+k1*delta_rv,zr);
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
                num1=abs(numerator)^2;
        num2=denominator;
        EDOF=(abs(numerator)^2)/denominator;

end