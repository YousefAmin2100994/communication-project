T=100
ts=0.001;
N=ceil(T/ts);
fs=1/ts;
df=1/T;

t=-(N/2)*ts:ts:((N-1)/2)*ts;
v=size(t);
x=zeros(size(t));
x(t>=-2&t<-1)=t(t>=-2&t<-1)+2;
x(t>=-1&t<1)=1;
x(t>=1&t<2)=-t(t>=1&t<2)+2;
figure(1)
plot(t,x);
xlabel('time(sec)');ylabel('x(t)');grid on;


%%xlabel('time(s)');ylabel('x(t)');grid on;
%%v=3*sinc(2*t).*sinc(3*t);

%%plot(t,v);


%%%%%%%%%%%% spectrum in freq. domain

if(rem(N,2)==0)
f=-0.5*fs:df:0.5*fs-df;
else
f=-(0.5*fs-0.5*df):df:0.5*fs-0.5*df
end

X=fftshift(fft(x))*0.001;
figure(2)
plot(f,abs(X));
xlabel('freq(HZ)');ylabel('X(f))');grid on;
%%%%%%%%%%%%%%
powersigt = sum(abs(x).^2)*0.001;
powersigf = sum(abs(X).^2)/T;




%for calculating bw that contain 95% of ths signal power


index=find(f==0);
power_acc = 0;
for c_index = index:length(f)
    power_acc = (abs(X(c_index)).^2)/T + power_acc;
    if (power_acc >= 0.95*0.5*powersigf);
        BW = abs(f(c_index));
        break
    end
end
H1=abs(f)<1;
figure(3)
plot(f,abs(H1));
XO=real(ifft(ifftshift(H1.*X)/ts));
figure(4)

plot(t,XO,'color','r',t,x,'color','b');
legend('after filter','before filter')
xlabel('time(s)'); ylabel('real massege after filter');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=zeros(size(t))
y(t>0&t<6)=cos(2*pi*t(t>0&t<6));
figure(5)
plot(t,y);title('massege before transmission');
Y_freq=fftshift(fft(y)*ts)
figure(6)
plot(f,abs(Y_freq))
index_2=find(f==0);
power_acc_2 = 0;
powersigt_2 = sum(abs(y).^2)*0.001;
powersigf_2 = sum(abs(Y_freq).^2)/T;
for c_index_2 = index_2:length(f)
    power_acc_2 = (abs(Y_freq(c_index_2)).^2)/T + power_acc_2;
    if (power_acc_2 >= 0.95*0.5*powersigf_2);
        BW_2 = abs(f(c_index_2));
        break
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Frequency division multiplixing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%TX DSB_SC
c1=cos(2*pi*20*t);
C=fftshift(fft(c1))/N;
figure(7);
plot(f,C);

figure(8)
plot(t,c1)
xlabel('time(s)');ylabel('carrier signal')
s1=XO.*c1;
figure(9)
plot(t,s1);
xlabel('time(s)');ylabel('modulated signal');
%%%%%%%%%%%%%%%%%%%%%% lets show its spectrum in freq domain
S1=fftshift(fft(s1))*ts;
figure(10);
plot(f,S1);
xlabel('freq(HZ)');ylabel('modulated signal');
%%%%%%%%%%%%%%%%%%%%bandpass filter at fc

H2=(abs(f) >= 17)&(abs(f) <= 23);
figure(11);
plot(f,H2);title('BANDPASS FLITER @fc1');
S1_DASH=H2.*S1;
figure(12);
plot(f,S1_DASH);
xlabel('freq(HZ)');ylabel('modulated signal');

%%%%%%%%%%%%%%%% 2nd signal m(t)
c2=cos(2*pi*24.5*t);
C2=fftshift(fft(c1))/N;
figure(13)
plot(t,c2)
figure(14)
plot(f,C2)
s2=y.*c2;
figure(15);
plot(t,s2);
xlabel('time(s)');ylabel('modulated signal');
S2=fftshift(fft(s2))*ts;
figure(16);
plot(f,S2);
xlabel('freq(HZ)');ylabel('modulated signal');
%%%%%%%%%%%%%%%%%%%%%%%%

H3=(abs(f) >= 23)&(abs(f) <= 24.5);
figure(17);
plot(f,H3);title('BANDPASS FLITER @fc2');
S2_DASH=H3.*S2;
figure(18);
plot(f,S2_DASH);
xlabel('freq(HZ)');ylabel('modulated signal');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=S1_DASH+S2_DASH;
figure(19);
plot(f,S);
H4=abs(f)<24.5;
st=real(ifft(ifftshift(H4.*S)))
figure(20);
plot(t,st);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DEMODULATION PROCESS
V1=S.*H2;
v1=real(ifft(ifftshift(V1)));
v1n=v1.*c1;
HL=abs(f)<1.5;
Xo_recieved=HL.*(fftshift(fft(v1n)));
xo_recieved=real(ifft(ifftshift(Xo_recieved))/ts);
figure(21);
plot(t,xo_recieved); title('massege recieved X(t)');

V2=S.*H3;
v2=real(ifft(ifftshift(V2)));
v2n=v2.*c2;
HL=abs(f)<1.5;
MO_recieved=HL.*(fftshift(fft(v2n)));
mo_recieved=real(ifft(ifftshift(MO_recieved))/ts);
figure(22);
plot(t,mo_recieved); title('massege recieved m(t)');
figure(23)
plot(t,XO);title('massege before transmission');















