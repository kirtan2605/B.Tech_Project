K=[0:0.1:30];
n=length(K);
for i=1:n
numg=[1 0]; deng=[1 10 10 5*K(i)];
sys = tf(numg,deng);
p(:,i)=pole(sys);
end
plot(real(p),imag(p),’.’), grid
xlabel(’real’);
ylabel(’imaginary’);
text(0,3.1263,’K=20’) ;
text(0.0,-3.1263,’K=20’);
