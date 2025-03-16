function [EDOF1]=sULAclosed(Nt,Nr,Lt,Lr,fc,c,xr,yr,zr,u,l,s,f,L,Ts,Channel,laisi_k)
        M=Nt;
        N=Nr;
        intervalt=Lt/M;
        intervalr=Lr/N;
        denominator=0;
         numerator=0; 
                 distances=sqrt(xr^2+yr^2+zr^2);
        [channelGaindB,~] = functionChannelgain(distances);
        channelGain=sqrt(db2pow(channelGaindB));
         for m=1:M
            for n =1:N
                 R_t_m=location(0,-Lt/2+(m-1)*intervalt,0);
               R_r_n=location(xr,yr-Lr/2+(n-1)*intervalr,zr);
               [w,~]=doppler(R_t_m,R_r_n,u,f(l));
               w=w*Ts;
               
               if l~=s
                   Channel(n,m)=Channel(n,m)*(-1)^(l-s)*w/sqrt(2)/(l-s);
               end

            end
         end
               R_connection=Channel'*Channel;
               trR=trace(R_connection);
               trR2=trace(R_connection^2);
        for m1=1:M
            for m2 =1:M
                num=0;
               for n=1:N
               R_t_m1=location(0,-Lt/2+(m1-1)*intervalt,0);
               R_t_m2=location(0,-Lt/2+(m2-1)*intervalt,0);
               R_r_n=location(xr,yr-Lr/2+(n-1)*intervalr,zr);
               d=sqrt(xr^2+(zr)^2);
               fd1=doppler(R_t_m1,R_r_n,u,f(l))*Ts;
                 k0=2*pi*f(l)/c;
               fd2=doppler(R_t_m2,R_r_n,u,f(l))*Ts;
               fd=(fd1+fd2)/2;
               Dnm2=sqrt(d^2+(intervalt*(m2-1)-yr-(n-1)*intervalr)^2);
               Dnm1=sqrt(d^2+(intervalt*(m1-1)-yr-(n-1)*intervalr)^2);
               if (l-s+fd)~=0
               num=num+channelGain*channelGain*(distances/Dnm1)*(distances/Dnm2)*exp(-1j*k0*Dnm2+1j*k0*Dnm1)*sin(pi*(l-s+fd))^2/(L^2*sin(pi*(l-s+fd)/L)^2);
%              num=num+exp(-1j*k0*Dnm2+1j*k0*Dnm1)*sin(pi*(l-s+fd1))*sin(pi*(l-s+fd2))/(L*sin(pi*(l-s+fd1)/L))/(L*sin(pi*(l-s+fd2)/L));

               else
               num=num+channelGain*channelGain*(distances/Dnm1)*(distances/Dnm2)*exp(-1j*k0*Dnm2+1j*k0*Dnm1);
               end
               %               num=num+exp(-1j*k0*Dnm2+1j*k0*Dnm1)*conj(DFO(l,s,fd1,L))*DFO(l,s,fd2,L);
               if m1==1
                   if (l-s+fd2)~=0
               numerator=numerator+channelGain^2*(distances/Dnm2)^2*sin(pi*(l-s+fd2))^2/(L^2*sin(pi*(l-s+fd2)/L)^2);
                   else
                numerator=numerator+channelGain^2*(distances/Dnm2)^2 ;      
%                numerator=numerator+1/(d^2+((zr+(n-1)*intervalr)-(m2-1)*intervalt)^2)*DFO(l,s,fd2,L)'*DFO(l,s,fd2,L);
                   end
                   end
               end
               num=abs(num)^2;
               denominator=denominator+num;
            end
        end
        EDOF1=(sqrt(laisi_k/(1+laisi_k))^2*abs(numerator)+trR)^2/(sqrt(laisi_k/(1+laisi_k))^4*denominator+trR2);

end