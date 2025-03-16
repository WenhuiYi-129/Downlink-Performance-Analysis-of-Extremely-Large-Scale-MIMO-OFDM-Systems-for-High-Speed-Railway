function ending=DFO(l,s,fd,L)

if (l-s+fd)~=0
ending=sin(pi*(l-s+fd))*exp(1j*pi*(1-1/L)*(l-s+fd))/(L*sin(pi*(l-s+fd)/L));
else
    ending=1;
end

end