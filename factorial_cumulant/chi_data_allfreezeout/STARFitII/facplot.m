clear;
c1=load("./chi1.dat");
c2=load("./chi2.dat");
c3=load("./chi3.dat");
c4=load("./chi4.dat");
sqrts=load("./sqrts.dat");
c1(11)=0.030;
c1(12)=0.033;
c1(13)=0.035;
c1(14)=0.038;
c1(15)=0.040;
c1(16)=0.042;
c1(17)=0.045;
c1(18)=0.048;
c1(19)=0.052;
c1(20)=0.055;
c2(11)=0.0950;
c2(12)=0.0954;
c2(13)=0.0959;
c2(14)=0.0964;
c2(15)=0.0970;
c2(16)=0.0977;
c2(17)=0.0984;
c2(18)=0.0992;
c2(19)=0.1000;
c2(20)=0.1009;
c1(1:40)=smooth(c1(1:40),20);
c2(1:40)=smooth(c2(1:40),20);
c3(1:40)=smooth(c3(1:40),30);
c4(1:40)=smooth(c4(1:40),30);
%plot(c2(1:40))
k1=c1;
k2=c2-c1;
k3=2*c1-3*c2+c3;
k4=-6*c1+11*c2-6*c3+c4;
Rk21=k2./k1;
Rk31=k3./k1;
Rk41=k4./k1;
r21=c2./c1;
R21HRG=load("./R21HRGII.dat")';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%semilogx(sqrts,c3./c1,'DisplayName','fRG-GCE');
%axis([1 200 0.1 3]);

%figure;
%semilogx(sqrts,c4./c2,'DisplayName','fRG-GCE');
%axis([1 200 0.1 3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R21star=load("../STAR_data/R21.dat");
R31star=load("../STAR_data/R31.dat");
R41star=load("../STAR_data/R41.dat");
R31starup=load("../STAR_data/R31_up.dat");
R41starup=load("../STAR_data/R41_up.dat");
x=R21star(:,1);
r21cen=R21star(:,2);
r31cen=R31star(:,2);
r41cen=R41star(:,2);
R31starerr=R31starup-R31star(:,2);
R41starerr=R41starup-R41star(:,2);

figure;
%semilogx(sqrts,r21-R21HRG);
%legend('$R^{fRG}_{21}-R^{HRG}_{21}$', 'Interpreter', 'latex');
%hold on;
semilogx(sqrts,r21,'DisplayName','fRG');
hold on;
axis([1 250 -0.2 10]);
semilogx(sqrts,R21HRG,'DisplayName','HRG');
semilogx(x,r21cen,'o','DisplayName','STAR-proton');
legend show
xlabel('$\sqrt{s_{NN}}$', 'Interpreter', 'latex')
ylabel('$R^f_{21}$', 'Interpreter', 'latex')
hold off;

%figure;
%semilogx(sqrts,Rk31,'DisplayName','fRG-GCE');
%axis([0 200 -0.2 0.8]);
%hold on;
%errorbar(x,r31cen,R31starerr,'o','DisplayName','STAR-proton')
%legend show
%xlabel('$\sqrt{s_{NN}}$', 'Interpreter', 'latex')
%ylabel('$R^f_{31}$', 'Interpreter', 'latex')
%hold off;

%figure;
%semilogx(sqrts,Rk41,'DisplayName','fRG-GCE');
%axis([0 200 -0.5 3]);
%hold on;
%errorbar(x,r41cen,R41starerr,'o','DisplayName','STAR-proton')
%legend show
%xlabel('$\sqrt{s_{NN}}$', 'Interpreter', 'latex')
%ylabel('$R^f_{41}$', 'Interpreter', 'latex')
%hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%