function[EDOF1]= sUPAclosed(Nth,Ntv,Nrh,Nrv,L_tH,L_tV,L_rH,L_rV,fc,c,xr,yr,zr,u,l,s,f,L,Ts,Channel,laisi_k)
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
        denominator=0;
         numerator=0; 
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
% trR=0;
% trR2=0;
        for m1=1:M
             i1=mod(m1-1,TX_number_h);
             j1=floor((m1-1)/TX_number_h);
            for m2 =1:M
                i2=mod(m2-1,TX_number_h);
                j2=floor((m2-1)/TX_number_h);
                num=0;
               for n=1:N
                t1=mod(n-1,RX_number_h);
                k1=floor((n-1)/RX_number_h);
               R_t_m1=location(-L_tH/2+i1*delta_th,-L_tV/2+j1*delta_tv,0);
               R_t_m2=location(-L_tH/2+i2*delta_th,-L_tV/2+j2*delta_tv,0);
               R_r_n=location(xr-L_rH/2+t1*delta_rh,yr-L_rV/2+k1*delta_rv,zr);
               d=sqrt(xr^2+(zr)^2);
               [~,fd1]=doppler(R_t_m1,R_r_n,u,f(l));
               
               fd1=fd1*Ts;
               k0=2*pi*f(l)/c;
               [~,fd2]=doppler(R_t_m2,R_r_n,u,f(l));
             
               fd2=fd2*Ts;
               fd=(fd1+fd2)/2;
               Dnm1=sqrt(zr^2+((-L_tH/2+i1*delta_th)-(xr-L_rH/2+t1*delta_rh))^2+((-L_tV/2+j1*delta_tv)-(yr-L_rV/2+k1*delta_rv))^2);
               Dnm2=sqrt(zr^2+((-L_tH/2+i2*delta_th)-(xr-L_rH/2+t1*delta_rh))^2+((-L_tV/2+j2*delta_tv)-(yr-L_rV/2+k1*delta_rv))^2);
               if (l-s+fd)~=0
               num=num+channelGain*channelGain*(distances/Dnm1)*(distances/Dnm2)*exp(-1j*k0*Dnm2+1j*k0*Dnm1)*sin(pi*(l-s+fd))^2/(L^2*sin(pi*(l-s+fd)/L)^2);
               else
                   num=num+channelGain*channelGain*(distances/Dnm1)*(distances/Dnm2)*exp(-1j*k0*Dnm2+1j*k0*Dnm1);
               end
               if m1==1
                  if (l-s+fd2)~=0
%                    numerator=numerator+(Dnm2/distances)^2*1/(xr^2+(-L_rH/2+yr+t1*delta_rh+L_tH/2-i2*delta_th)^2+(-L_rV/2+zr+k1*delta_rv+L_tV/2-j2*delta_tv)^2)*sin(pi*(l-s+fd2))^2/(L^2*sin(pi*(l-s+fd2)/L)^2);
                  numerator=numerator+channelGain^2*(distances/Dnm2)^2*sin(pi*(l-s+fd2))^2/(L^2*sin(pi*(l-s+fd2)/L)^2);
                  else
%                     numerator=numerator+(Dnm2/distances)^2*1/(xr^2+(-L_rH/2+yr+t1*delta_rh+L_tH/2-i2*delta_th)^2+(-L_rV/2+zr+k1*delta_rv+L_tV/2-j2*delta_tv)^2);
               numerator=numerator+channelGain^2*(distances/Dnm2)^2;

                  end
               end
               end
                num=abs(num)^2;
               denominator=denominator+num;
            end
        end
        num1=abs(numerator)^2;
        num2=denominator;
        EDOF1=(sqrt(laisi_k/(1+laisi_k))^2*abs(numerator)+trR)^2/(sqrt(laisi_k/(1+laisi_k))^4*denominator+trR2);


end