e=-40:2:40;
L=5;
for i=1:length(e)
T(i)=sin(pi*(1+e(i)))/L/sin(pi*(1+e(i))/L);
end
plot(e,T);


