%___________________________________________________________________%
%  FS-MPC predictive control of matrix converters                   %
%                                                                   %
%  Developed in MATLAB R2016a                                       %
%                                                                   %
%                         DJEOUADI MOHAMMED                         %
%                                                                   %
%         e-Mail: mohammeddjeouadi3@gmail.com                       %
%                                                                   %
%                                                                   %
%                                                                   %
%___________________________________________________________________%

clear all
clc
%sig; 
f=50;w=2*pi*f;k=1;Vm=690*sqrt(2);h=0.001;tf=1;n=1000;Ts=50e-6;
Rr=54.44e-3;Rs=4.45e-3;Lf=0.75e-3;Rf=75e-3;Cf=38e-3;Ls=134e-6;Lr=1.6e-3;p=2;M=3.33;sig=33;flis=55;Wr=w-p;Us=690;Im=20;
Psref=12e3;Qsref=0;
iqr(k)=0;idr(k)=0;
Vdr(1)=0;Vqr(1)=0;
Ac=[0 1/Cf;-1/Lf -Rf/Lf];Bc=[0 -1/Cf;1/Lf 0];
%Identity matrix 
I=eye(2);
Aq=exp(Ac*Ts); Bq=(Ac)'*(Aq-I)*Bc;
jj=1;
% Valid switching states of the Matrix Converter
S(:,:,1) = [1 0 0; 1 0 0; 1 0 0];
S(:,:,2) = [0 1 0; 0 1 0; 0 1 0];
S(:,:,3) = [0 0 1; 0 0 1; 0 0 1];
S(:,:,4) = [1 0 0; 0 0 1; 0 0 1];
S(:,:,5) = [0 1 0; 0 0 1; 0 0 1];
S(:,:,6) = [0 1 0; 1 0 0; 1 0 0];
S(:,:,7) = [0 0 1; 1 0 0; 1 0 0];
S(:,:,8) = [0 0 1; 0 1 0; 0 1 0];
S(:,:,9) = [1 0 0; 0 1 0; 0 1 0];
S(:,:,10) = [0 0 1; 1 0 0; 0 0 1];
S(:,:,11) = [0 0 1; 0 1 0; 0 0 1];
S(:,:,12) = [1 0 0; 0 1 0; 1 0 0];
S(:,:,13) = [1 0 0; 0 0 1; 1 0 0];
S(:,:,14) = [0 1 0; 0 0 1; 0 1 0];
S(:,:,15) = [0 1 0; 1 0 0; 0 1 0];
S(:,:,16) = [0 0 1; 0 0 1; 1 0 0];
S(:,:,17) = [0 0 1; 0 0 1; 0 1 0];
S(:,:,18) = [1 0 0; 1 0 0; 0 1 0];
S(:,:,19) = [1 0 0; 1 0 0; 0 0 1];
S(:,:,20) = [0 1 0; 0 1 0; 0 0 1];
S(:,:,21) = [0 1 0; 0 1 0; 1 0 0];
S(:,:,22) = [1 0 0; 0 1 0; 0 0 1];
S(:,:,23) = [1 0 0; 0 0 1; 0 1 0];
S(:,:,24) = [0 1 0; 1 0 0; 0 0 1];
S(:,:,25) = [0 1 0; 0 0 1; 1 0 0];
S(:,:,26) = [0 0 1; 1 0 0; 0 1 0];
S(:,:,27) = [0 0 1; 0 1 0; 1 0 0];
for k=1:n
    t=k*Ts;
Vs1(k)=Vm*sin(2*pi*f*t);
Vs2(k)=Vm*sin(2*pi*f*t-(2*pi/3));
Vs3(k)=Vm*sin(2*pi*f*t+(2*pi/3));
Vs=[Vs1;Vs2;Vs3];
Vr=(2/3)*[cos(w*t) cos(w*t-(2*pi/3)) cos(w*t+(2*pi/3));-sin(w*t) -sin(w*t-(2*pi/3)) -sin(w*t+(2*pi/3));0.5 0.5 0.5];
Vdq=Vr*Vs;
% Vd(k)=Vdq(1);
% Vq(k)=Vdq(2);
Vd(k)=(2/3)*(Vs1(k)*cos(w*t)+Vs2(k)*cos(w*t-(2*pi/3))+Vs3(k)*cos(w*t+(2*pi/3))); 
Vq(k)=(-2/3)*(Vs1(k)*sin(w*t)+Vs2(k)*sin(w*t-(2*pi/3))+Vs3(k)*sin(w*t+(2*pi/3)));

iqr(k+1)=(Ts/(sig*Lr))*(Vqr(k)-Rr*iqr(k)-Wr*sig*Lr*idr(k)-Wr*(M/Ls)*flis)+iqr(k);
idr(k+1)=(Ts/(sig*Lr))*(Vdr(k)-Rr*idr(k)+Wr*sig*Lr*iqr(k))+idr(k);

    iar(k+1)=idr(k+1)*cos(w*t)-iqr(k+1)*sin(w*t);
    ibr(k+1)=idr(k+1)*cos(w*t-(2*pi/3))-iqr(k+1)*sin(w*t-(2*pi/3));
    icr(k+1)=idr(k+1)*cos(w*t+(2*pi/3))-iqr(k+1)*sin(w*t+(2*pi/3));

Ps(k)=-(M/Ls)*Us*iqr(k+1);
Cem(k)=-p*(M/Ls)*flis*iqr(k+1);
Qs(k)=Us*flis/Ls-Us*(M/Ls)*idr(k+1);


% k=1:n;
% t=k*Ts;
% 
% plot(t,Vd(k),t,Vq(k))
% 
% plot(t,Ps(k),t,Cem(k),t,Qs(k))
   
      Va(k)=0;Vb(k)=0;Vc(k)=0;
    iqr(k+1)=0;idr(k+1)=0;iaf(k)=0;ibf(k)=0;icf(k)=0;


    iaf(k+1)=Aq(2,1)*Va(k)+Aq(2,2)*iaf(k)+Bq(2,1)*Vs1(k)+Bq(2,2)*iar(k+1);
    ibf(k+1)=Aq(2,1)*Vb(k)+Aq(2,2)*ibf(k)+Bq(2,1)*Vs2(k)+Bq(2,2)*ibr(k+1);
    icf(k+1)=Aq(2,1)*Vc(k)+Aq(2,2)*icf(k)+Bq(2,1)*Vs3(k)+Bq(2,2)*icr(k+1);
    
    Vaf(k+1)=Aq(1,1)*Va(k)+Aq(1,2)*iaf(k+1)+Bq(1,1)*Vs1(k)+Bq(1,2)*iar(k+1);
    Vbf(k+1)=Aq(1,1)*Vb(k)+Aq(1,2)*ibf(k+1)+Bq(1,1)*Vs2(k)+Bq(1,2)*ibr(k+1);
    Vcf(k+1)=Aq(1,1)*Vc(k)+Aq(1,2)*icf(k+1)+Bq(1,1)*Vs3(k)+Bq(1,2)*icr(k+1);
    

    %voltage predictions of network
    Vdf(k+1)=(2/3)*(Vaf(k+1)*cos(w*t)+Vbf(k+1)*cos(w*t-(2*pi/3))+Vcf(k+1)*cos(w*t+(2*pi/3))); 
    Vqf(k+1)=(-2/3)*(Vaf(k+1)*sin(w*t)+Vbf(k+1)*sin(w*t-(2*pi/3))+Vcf(k+1)*sin(w*t+(2*pi/3)));
    %current predictions of network
    idf(k+1)=(2/3)*(iaf(k+1)*cos(w*t)+ibf(k+1)*cos(w*t-(2*pi/3))+icf(k+1)*cos(w*t+(2*pi/3))); 
    iqf(k+1)=(-2/3)*(iaf(k+1)*sin(w*t)+ibf(k+1)*sin(w*t-(2*pi/3))+icf(k+1)*sin(w*t+(2*pi/3)));
    %The reactive power that circulates between the network and the rotor 
    Qf(k+1)=Vqf(k+1)*idf(k+1)-Vdf(k+1)*iqf(k+1);
    % input current vector is given by the switches state
    % with the load voltage
%     iA(jj)=S(1,1,jj)*ia(k)+S(1,2,jj)*ib(k)+S(1,3,jj)*ic(k);
%     iB(jj)=S(2,1,jj)*ia(k)+S(2,2,jj)*ib(k)+S(2,3,jj)*ic(k);
%     iC(jj)=S(3,1,jj)*ia(k)+S(3,2,jj)*ib(k)+S(3,3,jj)*ic(k);
    
    
    VA(k+1)=S(1,1,jj)*Vaf(k+1)+S(1,2,jj)*Vbf(k+1)+S(1,3,jj)*Vcf(k+1);
    VB(k+1)=S(2,1,jj)*Vaf(k+1)+S(2,2,jj)*Vbf(k+1)+S(2,3,jj)*Vcf(k+1);
    VC(k+1)=S(3,1,jj)*Vaf(k+1)+S(3,2,jj)*Vbf(k+1)+S(3,3,jj)*Vcf(k+1);
    
    Vdp(k+1)=(2/3)*(VA(k+1)*cos(w*t)+VB(k+1)*cos(w*t-(2*pi/3))+VC(k+1)*cos(w*t+(2*pi/3))); 
    Vqp(k+1)=(-2/3)*(VA(k+1)*sin(w*t)+VB(k+1)*sin(w*t-(2*pi/3))+VC(k+1)*sin(w*t+(2*pi/3)));
    %Predictions of rotor currents
    iqrp(k+1)=(Ts/(sig*Lr))*(Vqp(k+1)-Rr*iqr(k+1)-Wr*sig*Lr*idr(k+1)-Wr*(M/Ls)*flis)+iqr(k+1);
    idrp(k+1)=(Ts/(sig*Lr))*(Vdp(k+1)-Rr*idr(k+1)+Wr*sig*Lr*iqr(k+1))+idr(k+1);
    Psp(k+1)=-(M/Ls)*Us*iqrp(k+1);
    Qsp(k+1) =Us*flis/Ls-Us*(M/Ls)*idrp(k+1);
   
     Va(k)=0;Vb(k)=0;Vc(k)=0;
    iqr(jj)=0;idr(jj)=0;iaf(k)=0;ibf(k)=0;icf(k)=0;

% Initialization of the optimal value of the cost function
gobt=inf;
% Calculation of predictions for the 27 switching states
for jj=1:27
%     t=jj*Ts;
%     Va(k)=0;Vb(k)=0;Vc(k)=0;
%     iqr(jj)=0;idr(jj)=0;iaf(k)=0;ibf(k)=0;icf(k)=0;
 
    
    %Inverse park transformation(dq0===>abc)

    %Va=Via..
    iaf(jj)=Aq(2,1)*Va(k)+Aq(2,2)*iaf(k)+Bq(2,1)*Vs1(k)+Bq(2,2)*iar(k+1);
    ibf(jj)=Aq(2,1)*Vb(k)+Aq(2,2)*ibf(k)+Bq(2,1)*Vs2(k)+Bq(2,2)*ibr(k+1);
    icf(jj)=Aq(2,1)*Vc(k)+Aq(2,2)*icf(k)+Bq(2,1)*Vs3(k)+Bq(2,2)*icr(k+1);
    
    Vaf(k+1)=Aq(1,1)*Va(k)+Aq(1,2)*iaf(jj)+Bq(1,1)*Vs1(k)+Bq(1,2)*iar(k+1);
    Vbf(k+1)=Aq(1,1)*Vb(k)+Aq(1,2)*ibf(jj)+Bq(1,1)*Vs2(k)+Bq(1,2)*ibr(k+1);
    Vcf(k+1)=Aq(1,1)*Vc(k)+Aq(1,2)*icf(jj)+Bq(1,1)*Vs3(k)+Bq(1,2)*icr(k+1);
    

    %voltage predictions of network
    Vdf(jj)=(2/3)*(Vaf(k+1)*cos(w*t)+Vbf(k+1)*cos(w*t-(2*pi/3))+Vcf(k+1)*cos(w*t+(2*pi/3))); 
    Vqf(jj)=(-2/3)*(Vaf(k+1)*sin(w*t)+Vbf(k+1)*sin(w*t-(2*pi/3))+Vcf(k+1)*sin(w*t+(2*pi/3)));
    %current predictions of network
    idf(jj)=(2/3)*(iaf(jj)*cos(w*t)+ibf(jj)*cos(w*t-(2*pi/3))+icf(jj)*cos(w*t+(2*pi/3))); 
    iqf(jj)=(-2/3)*(iaf(jj)*sin(w*t)+ibf(jj)*sin(w*t-(2*pi/3))+icf(jj)*sin(w*t+(2*pi/3)));
    %The reactive power that circulates between the network and the rotor 
    Qf(jj)=Vqf(jj)*idf(jj)-Vdf(jj)*iqf(jj);
    % input current vector is given by the switches state
    % with the load voltage
%     iA(jj)=S(1,1,jj)*ia(k)+S(1,2,jj)*ib(k)+S(1,3,jj)*ic(k);
%     iB(jj)=S(2,1,jj)*ia(k)+S(2,2,jj)*ib(k)+S(2,3,jj)*ic(k);
%     iC(jj)=S(3,1,jj)*ia(k)+S(3,2,jj)*ib(k)+S(3,3,jj)*ic(k);
    
    
    VA(jj)=S(1,1,jj)*Vaf(k+1)+S(1,2,jj)*Vbf(k+1)+S(1,3,jj)*Vcf(k+1);
    VB(jj)=S(2,1,jj)*Vaf(k+1)+S(2,2,jj)*Vbf(k+1)+S(2,3,jj)*Vcf(k+1);
    VC(jj)=S(3,1,jj)*Vaf(k+1)+S(3,2,jj)*Vbf(k+1)+S(3,3,jj)*Vcf(k+1);
    
    Vdp(jj)=(2/3)*(VA(jj)*cos(w*t)+VB(jj)*cos(w*t-(2*pi/3))+VC(jj)*cos(w*t+(2*pi/3))); 
    Vqp(jj)=(-2/3)*(VA(jj)*sin(w*t)+VB(jj)*sin(w*t-(2*pi/3))+VC(jj)*sin(w*t+(2*pi/3)));
    %Predictions of rotor currents
    iqrp(jj)=(Ts/(sig*Lr))*(Vqp(jj)-Rr*iqr(k+1)-Wr*sig*Lr*idr(k+1)-Wr*(M/Ls)*flis)+iqr(k+1);
    idrp(jj)=(Ts/(sig*Lr))*(Vdp(jj)-Rr*idr(k+1)+Wr*sig*Lr*iqr(k+1))+idr(k+1);
    Psp(jj)=-(M/Ls)*Us*iqrp(jj);
    Qsp(jj)=Us*flis/Ls-Us*(M/Ls)*idrp(jj);
    % cost function calculate
     g=abs(Psp(jj)-Psref)+abs(Qsp(jj)-Qsref);
    % optimization
    if g<gobt
        gobt=g;
        jobt=jj;
    end
    
end
 Vdr(k+1)=Vdp(jobt);
 Vqr(k+1)=Vqp(jobt);
end

k=1:n;
t=k*Ts;

% jj=1:27;
% t1=jj*Ts;


%I named (t1) because we already have (t) before ,
%if we use (t) in the same plot it can't show the graph because it's not the same lentgh 
% plot(t,Ps(k),t,Cem(k),t,Qs(k)),grid on,hold on
% plot(t1,Psp(jj),t1,Qsp(jj),t1,g)


plot(t,Psp(k+1),t,Qsp(k+1),t,g,'*'),grid on,hold on



