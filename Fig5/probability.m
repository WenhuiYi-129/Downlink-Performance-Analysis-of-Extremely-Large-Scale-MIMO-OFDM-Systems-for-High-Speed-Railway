function [ending]=probability(xr,zr,fc,Ts,dtransmit_position,itransmit1_position,itransmit2_position,Lth,Ltv,Lrh,Lrv,Mh,Mv,Nh,Nv,u,L,power,area,SINRth)
        
        num=0;
        Number=2*area;%%采样数Number of samples
        for i=1:Number
        yr1(i)=-area+i*2*area/Number;
        end
        for j=1:length(yr1)
            yr=yr1(j);
            receive_position=[xr yr zr];
            SINR1=SINR(fc,Ts,dtransmit_position,itransmit1_position,itransmit2_position,receive_position,Lth,Ltv,Lrh,Lrv,Mh,Mv,Nh,Nv,u,L,power);
            if SINR1>=SINRth
                num=num+1;
            end
        end        
        ending=num/Number;

end