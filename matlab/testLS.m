clear
a=0;
b=2;
N=100;
mu=10*N;
f=zeros(N,1);
B=zeros(N);
n=linspace(a,b,N+1);
h=(b-a)/N;
x=n(1:end-1)+h.*0.5;
L=zeros(N,3,2);
L(1,1,1)=1/(x(1)-x(2));
L(1,1,2)=-x(2)/(x(1)-x(2));
W=(2*(- x(1)^2 + x(1)*x(2) + x(1)*x(3) - x(2)^2 + x(2)*x(3) - x(3)^2));
L(1,2,1)=(x(2)-2*x(1)+x(3))/W;
L(1,2,2)=(-x(2)^2 + x(1)*x(2) - x(3)^2 + x(1)*x(3))/W;
L(2,1,1)=1/(x(2)-x(1));
L(2,1,2)=-x(1)/(x(2)-x(1));
L(2,2,1)=(x(1) - 2*x(2) + x(3))/W;
L(2,2,2)=(-x(1)^2 + x(2)*x(1) - x(3)^2 + x(2)*x(3))/W;
W=(2*(- x(2)^2 + x(2)*x(3) + x(2)*x(4) - x(3)^2 + x(3)*x(4) - x(4)^2));
L(2,3,1)=(x(3) - 2*x(2) + x(4))/W;
L(2,3,2)=(-x(3)^2 + x(2)*x(3) - x(4)^2 + x(2)*x(4))/W;
W=2*(x(1:end-4).^2+x(2:end-3).^2+x(3:end-2).^2-x(1:end-4).*x(2:end-3)-x(1:end-4).*x(3:end-2)-x(2:end-3).*x(3:end-2));
L(3:end-2,1,1)=(2*x(3:end-2)-x(2:end-3)-x(1:end-4))./W;
L(3:end-2,1,2)=(x(1:end-4).^2+x(2:end-3).^2-x(1:end-4).*x(3:end-2)-x(2:end-3).*x(3:end-2))./W;
W=-2*(x(2:end-3).^2+x(3:end-2).^2+x(4:end-1).^2-x(2:end-3).*x(3:end-2)-x(2:end-3).*x(4:end-1)-x(3:end-2).*x(4:end-1));
L(3:end-2,2,1)=(-2*x(3:end-2)+x(2:end-3)+x(4:end-1))./W;
L(3:end-2,2,2)=(-x(4:end-1).^2-x(2:end-3).^2+x(4:end-1).*x(3:end-2)+x(2:end-3).*x(3:end-2))./W;
W=-2*(x(3:end-2).^2+x(4:end-1).^2+x(5:end).^2-x(3:end-2).*x(4:end-1)-x(3:end-2).*x(5:end)-x(4:end-1).*x(5:end));
L(3:end-2,3,1)=(-2*x(3:end-2)+x(4:end-1)+x(5:end))./W;
L(3:end-2,3,2)=(-x(4:end-1).^2-x(5:end).^2+x(4:end-1).*x(3:end-2)+x(5:end).*x(3:end-2))./W;
W=-2*((x(end-3))^2+x(end-2)^2+x(end-1)^2-x(end-3)*x(end-2)-x(end-3)*x(end-1)-x(end-2)*x(end-1));
L(end-1,1,1)=(x(end-3) + x(end-2) - 2*x(end-1))/W;
L(end-1,1,2)=(- x(end-3)^2 + x(end-1)*x(end-3) - x(end-2)^2 + x(end-1)*x(end-2))/W;
W=-2*((x(end-2))^2+x(end-1)^2+x(end)^2-x(end-2)*x(end-1)-x(end-2)*x(end)-x(end-1)*x(end));
L(end-1,2,1)=(x(end-2) + x(end) - 2*x(end-1))/W;
L(end-1,2,2)=(- x(end-2)^2 + x(end-1)*x(end-2) - x(end)^2 + x(end-1)*x(end))/W;
L(end-1,3,1)=-1/(x(end)-x(end-1));
L(end-1,3,2)=x(end)/(x(end)-x(end-1));
L(end,1,1)=(x(end-2) + x(end-1) - 2*x(end))/W;
L(end,1,2)=(- x(end-2)^2 + x(end)*x(end-2) - x(end-1)^2 + x(end)*x(end-1))/W;
L(end,2,1)=1/(x(end)-x(end-1));
L(end,2,2)=-x(end-1)/(x(end)-x(end-1));
hold on
axis equal
axis([a-0.1 b+0.1 -1.5 1.5])
xh=linspace(n(1),n(2),10);
plot(xh,L(1,1,1)*xh+L(1,1,2),'b');
xh=linspace(n(2),n(3),10);
plot(xh,L(1,2,1)*xh+L(1,2,2),'b');
for i=2:N-1
    xh=linspace(n(i-1),n(i),10);
    plot(xh,L(i,1,1)*xh+L(i,1,2),'b');
    xh=linspace(n(i),n(i+1),10);
    plot(xh,L(i,2,1)*xh+L(i,2,2),'b');
    xh=linspace(n(i+1),n(i+2),10);
    plot(xh,L(i,3,1)*xh+L(i,3,2),'b');
end
xh=linspace(n(N-1),n(N),10);
plot(xh,L(N,1,1)*xh+L(N,1,2),'b');
xh=linspace(n(N),n(N+1),10);
plot(xh,L(N,2,1)*xh+L(N,2,2),'b');

y=sin(x);
xh=linspace(n(1),n(2),10);
yh=(L(1,1,1)*xh+L(1,1,2))*y(1)+(L(2,1,1)*xh+L(2,1,2))*y(2);
plot(xh,yh,'b');
plot(xh,sin(xh)-yh,'r');
xh=linspace(n(2),n(3),10);
yh=(L(1,2,1)*xh+L(1,2,2))*y(1);
yh=yh+(L(2,2,1)*xh+L(2,2,2))*y(2);
yh=yh+(L(3,1,1)*xh+L(3,1,2))*y(3);
plot(xh,yh,'b');
plot(xh,sin(xh)-yh,'r');
for i=3:N-2
    xh=linspace(n(i),n(i+1),10);
    yh=(L(i-1,3,1)*xh+L(i-1,3,2))*y(i-1);
    yh=yh+(L(i,2,1)*xh+L(i,2,2))*y(i);
    yh=yh+(L(i+1,1,1)*xh+L(i+1,1,2))*y(i+1);
    plot(xh,yh,'b');
    plot(xh,sin(xh)-yh,'r');
end
xh=linspace(n(end-2),n(end-1),10);
yh=(L(end-2,3,1)*xh+L(end-2,3,2))*y(end-2);
yh=yh+(L(end-1,2,1)*xh+L(end-1,2,2))*y(end-1);
yh=yh+(L(end,1,1)*xh+L(end,1,2))*y(end);
plot(xh,yh,'b');
plot(xh,sin(xh)-yh,'r');
xh=linspace(n(end-1),n(end),10);
yh=(L(end-1,3,1)*xh+L(end-1,3,2))*y(end-1);
yh=yh+(L(end,2,1)*xh+L(end,2,2))*y(end);
plot(xh,yh,'b');
plot(xh,sin(xh)-yh,'r');
%xh=linspace(a,b,N*9+1);
%plot(xh,sin(xh),'r');
%plot(x,y,'ro');

syms x;

%%
%right side
f(1)=int((L(1,1,1)*x+L(1,1,2))*cos(pi*x),n(1),n(2))+int((L(1,2,1)*x+L(1,2,2))*cos(pi*x),n(2),n(3));
for i=2:N-1
    f(i)=int((L(i,1,1)*x+L(i,1,2))*cos(pi*x),n(i-1),n(i))+int((L(i,2,1)*x+L(i,2,2))*cos(pi*x),n(i),n(i+1))+...
        int((L(i,3,1)*x+L(i,3,2))*cos(pi*x),n(i+1),n(i+2));
end
f(N)=int((L(N,1,1)*x+L(N,1,2))*cos(pi*x),n(N-1),n(N))+int((L(N,2,1)*x+L(N,2,2))*cos(pi*x),n(N),n(N+1));
%%
W1=((L(1,1,1)*n(2)+L(1,1,2))-(L(1,2,1)*n(2)+L(1,2,2)));
W2=L(1,2,1)*n(3)+L(1,2,2);
U1=0.5*(L(1,1,1)+L(1,2,1));
U2=0.5*L(1,2,1);
B(1,1)=int(L(1,1,1)*L(1,1,1),x,n(1),n(2))+int(L(1,2,1)*L(1,2,1),x,n(2),n(3))-...
    W1*U1-U1*W1-W2*U2- U2*W2+mu*W1*W1+mu*W2*W2;

W1=(L(2,1,1)*n(2)+L(2,1,2))-(L(2,2,1)*n(2)+L(2,2,2));
W2=(L(2,2,1)*n(3)+L(2,2,2))-(L(2,3,1)*n(3)+L(2,3,2));
W3=L(2,3,1)*n(4)+L(2,3,2);
U1=0.5*(L(2,1,1)+L(2,2,1));
U2=0.5*(L(2,2,1)+L(2,3,1));
U3=0.5*L(2,3,1);
B(2,2)=int(L(2,1,1)*L(2,1,1),x,n(1),n(2))+int(L(2,2,1)*L(2,2,1),x,n(2),n(3))+int(L(2,3,1)*L(2,3,1),x,n(3),n(4))-...
    W1*U1-U1*W1-W2*U2-U2*W2-W3*U3-U3*W3+...
    mu*W1*W1+mu*W2*W2+mu*W3*W3;
for i=3:N-2
   W1=-(L(i,1,1)*n(i-1)+L(i,1,2));
   W2=(L(i,1,1)*n(i)+L(i,1,2))-(L(i,2,1)*n(i)+L(i,2,2));
   W3=(L(i,2,1)*n(i+1)+L(i,2,2))-(L(i,3,1)*n(i+1)+L(i,3,2));
   W4=L(i,3,1)*n(i+2)+L(i,3,2);
   U1=0.5*L(i,1,1);
   U2=0.5*(L(i,1,1)+L(i,2,1));
   U3=0.5*(L(i,2,1)+L(i,3,1));
   U4=0.5*L(i,3,1);
   B(i,i)=int(L(i,1,1)*L(i,1,1),x,n(i-1),n(i))+int(L(i,2,1)*L(i,2,1),x,n(i),n(i+1))+int(L(i,3,1)*L(i,3,1),x,n(i+1),n(i+2))-...
    W1*U1-U1*W1-W2*U2-U2*W2-W3*U3-U3*W3-W4*U4-U4*W4+...
    mu*W1*W1+mu*W2*W2+mu*W3*W3+mu*W4*W4;
end
W1=-(L(N-1,1,1)*n(N-2)+L(N-1,1,2));
W2=(L(N-1,1,1)*n(N-1)+L(N-1,1,2))-(L(N-1,2,1)*n(N-1)+L(N-1,2,2));
W3=(L(N-1,2,1)*n(N)+L(N-1,2,2))-(L(N-1,3,1)*n(N)+L(N-1,3,2));
U1=0.5*L(N-1,1,1);
U2=0.5*(L(N-1,1,1)+L(N-1,2,1));
U3=0.5*(L(N-1,2,1)+L(N-1,3,1));
B(N-1,N-1)=int(L(N-1,1,1)*L(N-1,1,1),x,n(N-2),n(N-1))+int(L(N-1,2,1)*L(N-1,2,1),x,n(N-1),n(N))+int(L(N-1,3,1)*L(N-1,3,1),x,n(N),n(N+1))-...
    W1*U1-U1*W1-W2*U2-U2*W2-W3*U3-U3*W3+...
    mu*W1*W1+mu*W2*W2+mu*W3*W3;
W1=-(L(N,1,1)*n(N-1)+L(N,1,2));
W2=((L(N,1,1)*n(N)+L(N,1,2))-(L(N,2,1)*n(N)+L(N,2,2)));
U1=0.5*L(N,1,1);
U2=0.5*(L(N,1,1)+L(N,2,1));
B(N,N)=int(L(N,1,1)*L(N,1,1),x,n(N-1),n(N))+int(L(N,2,1)*L(N,2,1),x,n(N),n(N+1))-...
    W1*U1-U1*W1-W2*U2- U2*W2+mu*W1*W1+mu*W2*W2;
%%
%one point overlap
W1=L(1,2,1)*n(3)+L(1,2,2);
U1=0.5*L(1,2,1);
W12=-(L(4,1,1)*n(3)+L(4,1,2));
U12=0.5*L(4,1,1);
B(1,4)=-W1*U12-U1*W12+mu*W1*W12;
B(4,1)=B(1,4);
for i=2:N
    for j=2:N
        if(j-i==3)
            W1=L(i,3,1)*n(i+2)+L(i,3,2);
            U1=0.5*L(i,3,1);
            W12=-(L(j,1,1)*n(j-1)+L(j,1,2));
            U12=0.5*L(j,1,1);
            B(i,j)=-W1*U12-U1*W12+mu*W1*W12;
            B(j,i)=B(i,j);
        end
    end
end
%%
%two points overlap
W1=((L(1,1,1)*n(2)+L(1,1,2))-(L(1,2,1)*n(2)+L(1,2,2)));
W2=L(1,2,1)*n(3)+L(1,2,2);
U1=0.5*(L(1,1,1)+L(1,2,1));
U2=0.5*L(1,2,1);
W12=-(L(3,1,1)*n(2)+L(3,1,2));
W22=(L(3,1,1)*n(3)+L(3,1,2))-(L(3,2,1)*n(3)+L(3,2,2));
U12=0.5*L(3,1,1);
U22=0.5*(L(3,1,1)+L(3,2,1));
B(1,3)=int(L(1,2,1)*L(3,1,1),x,n(2),n(3))-...
    W1*U12-U1*W12-W2*U22-U2*W22+...
    mu*W1*W12+mu*W2*W22;
B(3,1)=B(1,3);
for i=2:N
    for j=2:N
        if(j-i==2)
           W1=((L(i,2,1)*n(i+1)+L(i,2,2))-(L(i,3,1)*n(i+1)+L(i,3,2)));
           W2=L(i,3,1)*n(i+2)+L(i,3,2);
           U1=0.5*(L(i,2,1)+L(i,3,1));
           U2=0.5*L(i,3,1);
           W12=-(L(j,1,1)*n(j-1)+L(j,1,2));
           W22=(L(j,1,1)*n(j)+L(j,1,2))-(L(j,2,1)*n(j)+L(j,2,2));
           U12=0.5*L(j,1,1);
           U22=0.5*(L(j,1,1)+L(j,2,1));
           B(i,j)=int(L(i,3,1)*L(j,1,1),x,n(i+1),n(i+2))-...
               W1*U12-U1*W12-W2*U22-U2*W22+mu*W1*W12+mu*W2*W22;
           B(j,i)=B(i,j);
        end
    end
end
%%
%three points overlap
W1=((L(1,1,1)*n(2)+L(1,1,2))-(L(1,2,1)*n(2)+L(1,2,2)));
W2=L(1,2,1)*n(3)+L(1,2,2);
U1=0.5*(L(1,1,1)+L(1,2,1));
U2=0.5*L(1,2,1);
W12=(L(2,1,1)*n(2)+L(2,1,2))-(L(2,2,1)*n(2)+L(2,2,2));
W22=(L(2,2,1)*n(3)+L(2,2,2))-(L(2,3,1)*n(3)+L(2,3,2));
U12=0.5*(L(2,1,1)+L(2,2,1));
U22=0.5*(L(2,2,1)+L(2,3,1));
B(1,2)=int(L(1,1,1)*L(2,1,1),x,n(1),n(2))+int(L(1,2,1)*L(2,2,1),x,n(2),n(3))-...
    W1*U12-U1*W12-W2*U22-U2*W22+...
    mu*W1*W12+mu*W2*W22;
B(2,1)=B(1,2);
for i=2:N-2
    for j=2:N-1
        if(j-i==1)
           W1=((L(i,1,1)*n(i)+L(i,1,2))-(L(i,2,1)*n(i)+L(i,2,2)));
           W2=((L(i,2,1)*n(i+1)+L(i,2,2))-(L(i,3,1)*n(i+1)+L(i,3,2)));
           W3=L(i,3,1)*n(i+2)+L(i,3,2);
           U1=0.5*(L(i,1,1)+L(i,2,1));
           U2=0.5*(L(i,2,1)+L(i,3,1));
           U3=0.5*L(i,3,1);
           W12=-(L(j,1,1)*n(j-1)+L(j,1,2));
           W22=(L(j,1,1)*n(j)+L(j,1,2))-(L(j,2,1)*n(j)+L(j,2,2));
           W32=(L(j,2,1)*n(j+1)+L(j,2,2))-(L(j,3,1)*n(j+1)+L(j,3,2));
           U12=0.5*L(j,1,1);
           U22=0.5*(L(j,1,1)+L(j,2,1));
           U32=0.5*(L(j,2,1)+L(j,3,1));
           B(i,j)=int(L(i,2,1)*L(j,1,1),x,n(i+1),n(i+2))+int(L(i,3,1)*L(j,2,1),x,n(i+2),n(i+3))-...
               W1*U12-U1*W12-W2*U22-U2*W22-W3*U32-U3*W32+...
               mu*W1*W12+mu*W2*W22+mu*W3*W32;
           B(j,i)=B(i,j);
        end
    end
end
 W1=((L(N-1,1,1)*n(N-1)+L(N-1,1,2))-(L(N-1,2,1)*n(N-1)+L(N-1,2,2)));
 W2=((L(N-1,2,1)*n(N)+L(N-1,2,2))-(L(N-1,3,1)*n(N)+L(N-1,3,2)));
 U1=0.5*(L(N-1,1,1)+L(N-1,2,1));
 U2=0.5*(L(N-1,2,1)+L(N-1,3,1));
 W12=-(L(N,1,1)*n(N-1)+L(N,1,2));
 W22=(L(N,1,1)*n(N)+L(N,1,2))-(L(N,2,1)*n(N)+L(N,2,2));
 U12=0.5*L(N,1,1);
 U22=0.5*(L(N,1,1)+L(N,2,1));
 B(N,N-1)=int(L(N-1,2,1)*L(N,1,1),x,n(N-1),n(N))+int(L(N-1,3,1)*L(N,2,1),x,n(N),n(N+1))-...
               W1*U12-U1*W12-W2*U22-U2*W22+...
               mu*W1*W12+mu*W2*W22;
 B(N-1,N)=B(N,N-1);
 %%
 %boundary 
  x=n(1:end-1)+h.*0.5;
  f2=pi^2*f;
 f2(1)=B(1,1)*cos(pi*x(1));
 f2(N)=B(N,N)*cos(pi*x(N));
 y(1)=cos(pi*x(1));
 y(N)=cos(pi*x(N));
 for i=2:4
     B(1,i)=0;
     f2(i)=f2(i)-B(i,1)*y(1)
     B(i,1)=0;
 end
 for i=N-3:N-1
     B(N,i)=0;
     f2(i)=f2(i)-B(i,N)*y(N)
     B(i,N)=0;
 end
 

 y1=cos(pi*x)'
 
 y=B\f2
 err=norm(y-y1)/norm(y1)