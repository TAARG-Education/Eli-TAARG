function [f1,c_D14,RatioD,xs,ys,D14,D34] = Universale2(Di,De,ch)
n=6;                      % grado del polinomio interpolante

% disp(['x = ',num2str(xh), ' | y = ',num2str(xv)]);

x1=([0.14545 0.25909 0.4     1       1.5545   1.7727   2]);
y1=([0.14545 0.21818 0.27273 0.16364 0.072727 0.045455 0.036364]);
xfit1=linspace(x1(1),x1(end),1000);
c=polyfit(x1,y1,n);
yfit02= polyval(c,xfit1);

y2=([0.36364 0.54545 0.61818 0.32727 0.12727 0.090909 0.072727]);
c=polyfit(x1,y2,n);
yfit04= polyval(c,xfit1);

y3=([0.85455 1.1273 1.1455 0.45455 0.18182 0.12727 0.1]);
c=polyfit(x1,y3,n);
yfit06= polyval(c,xfit1);

x4=([x1(1) 0.18182 x1(2) x1(3)  0.79091 x1(4)   1.0909  1.3818  x1(5)   x1(6)   1.8636  x1(7)]);
y4=([2.5   2.5636  2.4   1.8273 0.8     0.54545 0.45455 0.27273 0.21818 0.14909 0.12727 0.11273]);
xfit4=linspace(x4(1),x4(end),1000);
c=polyfit(x4,y4,11);
yfit08= polyval(c,xfit4);

x5=([0.3     0.31818 0.33182 0.35455 x1(3) 0.79091 x1(4)   1.0909  1.3818  x1(5)   x1(6)   1.8636  x1(7)]);
y5=([3.0000  2.8727  2.7091  2.4909  2.2   0.83636 0.56545 0.49091 0.28182 0.22727 0.16182 0.14545 0.12364]);
xfit5=linspace(x5(1),x5(end),1000);
c=polyfit(x5,y5,9);
yfit10= polyval(c,xfit5);

figure
plot(xfit1,yfit02,'k','LineWidth',1.5); hold on;
plot(xfit1,yfit04,'b','LineWidth',1.5);
plot(xfit1,yfit06,'r','LineWidth',1.5);
plot(xfit4,yfit08,'g','LineWidth',1.5);
plot(xfit5,yfit10,'m','LineWidth',1.5);
title('Velocità indotta dal vortice ad anello');
xlabel('c/D_{1/4}'); ylabel('f(c/D_{1/4},D_{3/4}/D_{1/4})');
% legend('D_{3/4}/D_{1/4} = 0.2','D_{3/4}/D_{1/4} = 0.4','D_{3/4}/D_{1/4} = 0.6','D_{3/4}/D_{1/4} = 0.8','D_{3/4}/D_{1/4} = 1.0');
text(x1(1)-0.1,yfit02(1) ,'0.2'); text(x1(1)-0.1,yfit04(1) ,'0.4'); text(x1(1)-0.1,yfit06(1),'0.6'); text(x1(1)-0.1,yfit08(1),'0.8'); text(0.5,2.7,'D_{3/4}/D_{1/4} = 1.0');
plot(x1,y1,'ko',x1,y2,'bo',x1,y3,'ro',x4,y4,'go',x5,y5,'mo','LineWidth',1);
axis ([0 2 0 3]);
grid minor;

for i=1:length(xfit1)
    yfit07(i)= ( yfit06(i)+yfit08(i) )/2;
end
plot(xfit1,yfit07,'k:','LineWidth',0.7); text(x1(1)-0.05,yfit07(1) ,'0.7');

for i=1:length(xfit1)
    yfit05(i)= ( yfit04(i)+yfit06(i) )/2;
end
plot(xfit1,yfit05,'k:','LineWidth',0.7); text(x1(1)-0.05,yfit05(1) ,'0.5');

for i=1:length(xfit1)
    yfit03(i)= ( yfit02(i)+yfit04(i) )/2;
end
plot(xfit1,yfit03,'k:','LineWidth',0.7); text(x1(1)-0.05,yfit03(1),'0.3');

for i=1:length(xfit1)
    yfit065(i)= ( yfit06(i)+yfit07(i) )/2;
end
plot(xfit1,yfit065,'k:','LineWidth',0.7); text(x1(1)-0.05,yfit065(1),'0.65');

for i=1:length(xfit1)
    yfit075(i)= ( yfit07(i)+yfit08(i) )/2;
end
plot(xfit1,yfit075,'k:','LineWidth',0.7); text(x1(1)-0.05,yfit075(1),'0.75');

for i=1:length(xfit1)
    yfit025(i)= ( yfit02(i)+yfit03(i) )/2;
end
plot(xfit1,yfit025,'k:','LineWidth',0.7); text(x1(1)-0.05,yfit025(1),'0.25');

for i=1:length(xfit1)
    yfit035(i)= ( yfit03(i)+yfit04(i) )/2;
end
plot(xfit1,yfit035,'k:','LineWidth',0.7); text(x1(1)-0.05,yfit035(1),'0.35');

for i=1:length(xfit1)
    yfit045(i)= ( yfit04(i)+yfit05(i) )/2;
end
plot(xfit1,yfit045,'k:','LineWidth',0.7); text(x1(1)-0.05,yfit045(1),'0.45');

for i=1:length(xfit1)
    yfit055(i)= ( yfit05(i)+yfit06(i) )/2;
end
plot(xfit1,yfit055,'k:','LineWidth',0.7); text(x1(1)-0.05,yfit055(1),'0.55');

for i=1:length(xfit1)
    yfit0625(i)= ( yfit06(i)+yfit065(i) )/2;
end
plot(xfit1,yfit0625,'k:','LineWidth',0.7); text(x1(1)-0.05,yfit0625(1),'0.625');

for i=1:length(xfit1)
    yfit0675(i)= ( yfit065(i)+yfit07(i) )/2;
end
plot(xfit1,yfit0675,'k:','LineWidth',0.7); text(x1(1)-0.05,yfit0675(1),'0.675');

for i=1:length(xfit1)
    yfit0725(i)= ( yfit07(i)+yfit075(i) )/2;
end
plot(xfit1,yfit0725,'k:','LineWidth',0.7); text(x1(1)-0.05,yfit0725(1),'0.725');

for i=1:length(xfit1)
    yfit0775(i)= ( yfit075(i)+yfit08(i) )/2;
end
plot(xfit1,yfit0775,'k:','LineWidth',0.7); text(x1(1)-0.05,yfit0775(1),'0.775');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t=linspace(0,1,100);
for i=1:length(t) 
xs(i)=0+t(i)*(ch);
ys(i)=Di/2+t(i)*(De/2-Di/2);
end
D14=2*ys(25);
D34=2*ys(75);

c_D14=ch/D14 ; % c/D_14
for i=1:length(xfit1)
    diff(i)=abs(c_D14-xfit1(i) );
    [m,j]=min(diff);
end
c_D14=xfit1(j);

CurvD=([0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.625 0.65 0.675 ,...
    0.7 0.725 0.75 0.775 0.8 1]);
YCurvD=([yfit02; yfit025; yfit03; yfit035; yfit04; yfit045; yfit05;,...
    yfit055; yfit06; yfit0625; yfit065; yfit0675; yfit07; yfit0725; yfit075;...
    yfit0775;yfit08;yfit10]);
RatioD=D34/D14;             %D_34/D_14
for i=1:length(CurvD)
    diffD(i)=abs(RatioD-CurvD(i) );
    [mi,k]=min(diffD);
end
RatioD=CurvD(k);
yy=ones(1,1000);
f1=YCurvD(k,j);
%Primo tentativo theta
mis_z=xs(end)-xs(75);
mis_ver=ys(75)-ys(end);
thetad=atand(mis_ver/mis_z);
theta=atan(mis_ver/mis_z);
end