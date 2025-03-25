%-------------------------------------------------------------------------
% Solving multi-term time fractional Convection Diffusion reaction Equation 
%-------------------------------------------------------------------------
%  (∂^β)V(x,t)−a(t)U_xx(x,t) + d(t)U_x(x,t) + κ(x,t)U(x,t) = f(x,t) 
%-------------------------------------------------------------------------
% using L1-scheme time (in graded mesh) & 4th order Compact Finite Differencce Method in space

clear all
tic
M =64; %grid points in space
N = 64; %grid points in time
X1 = 0;
X2 = 1;
T = 1.0;
dx = (X2-X1)/M;
%c=1;
b1=1;
b2=1;
b3=1;
b4=1;
% dt=T/N;

lamda1=0.4; %alpha lies in (0,1) % lamda2
lamda2=0.3;
lamda3=0.2;
lamda4=0.1;

%r=1; %(uniform case)
%r=(2-lamda1)/(2*lamda1);
%r=(2-lamda1)/(lamda1);
r=2*(2-lamda1)/(lamda1);

for j=1:N+1
    t(j)=T*((j-1)/N)^r;
    w(j)=j;   %time discretization 

    q(j)= (t(j))^4+t(j);         % q(t) coefficient of Uxx
     h(j)=0;             % P(t) coefficient of Ux
end

for i=1:M+1
    x(i)=X1+(i-1)*dx; %Space discretization
end

     for i=1:M+1
          for j=1:N+1
g(i,j)=exp(2*(t(j))-x(i)) ;     % g(x,t) coefficient of u
   g1(i,j)= 0;          % g1 is coefficients of g'(x,t)
    g2(i,j)= 0 ;         % g2 is coefficients of g''(x,t)
          end
end

%Initial condition
for i=1:M+1
    x(i)=X1+(i-1)*dx; %Space discretization
    phi(i)=0; %initial condition
    u(i,1)=phi(i);
end

for j=1:N
    dt(j)=t(j+1)-t(j); % non-uniform mesh
end

%Boundary conditions
for j=1:N+1
    leftbc(j)=0;
     rightbc(j)=0;
     
    u(1,j)=leftbc(j);    % boundary conditon 1
    u(M+1,j)=rightbc(j); % boundary conditon 2
end
P1=1/gamma(2-lamda1); 
P2=1/gamma(2-lamda2);
P3=1/gamma(2-lamda3);
P4=1/gamma(2-lamda4);


for j=1:N+1
    for  i=1:M+1 
        e(i,j)= t(j)^(lamda1)*(x(i)^2-x(i)^3);                                          % exact solution
        
        s(i,j)=b4*gamma(lamda1+1)/gamma(lamda1-lamda4+1)*t(j)^(lamda1-lamda4)*(x(i)^2-x(i)^3)+b3*gamma(lamda1+1)/gamma(lamda1-lamda3+1)*t(j)^(lamda1-lamda3)*(x(i)^2-x(i)^3)+b2*gamma(lamda1+1)/gamma(lamda1-lamda2+1)*t(j)^(lamda1-lamda2)*(x(i)^2-x(i)^3)+b1*(x(i)^2-x(i)^3)*lamda1*pi*csc(lamda1*pi)/gamma(1-lamda1)-((t(j))^((lamda1+4))+(t(j))^((lamda1+1)))*(2-6*x(i))+(exp(2*(t(j))-x(i)))*((t(j))^(lamda1))*((x(i))^2)*(1-x(i));            % R.H.S. function
    end
end

for k=1:N+1
    for i=2:M
        sum1(i,k)=0;
    end
end

for k=1:N+1
    for i=2:M
        sum2(i,k)=0;
    end
end

for k=1:N+1
    for i=2:M
        sum3(i,k)=0;
    end
end

for k=1:N+1
    for i=2:M
        sum4(i,k)=0;
    end
end


for j=1:N
    A=zeros(M+1,M+1);
    
    %First row elements by using the u(0,t)=t^(alpha)
    A(1,1)=1;
    
    %second row to mth row elements by using the gridpoints x_0 to x_m
    for i=2:M
            A(i,i-1) =(-q(j+1)/(dx)^2)+(1/12)*(g(i,j+1)-(((h(j+1))^2)/(q(j+1))^2))-(1/(2*dx))*h(j+1)-(((dx)/24)*(2*(g1(i,j+1))-((h(j+1)*g(i,j+1))/q(j+1))))-b1*P1*((dt(j))^(-lamda1))*((((-h(j+1)*dx)/q(j+1))-2)/24)-b2*P2*((dt(j))^(-lamda2))*((((-h(j+1)*dx)/q(j+1))-2)/24)-b3*P3*((dt(j))^(-lamda3))*((((-h(j+1)*dx)/q(j+1))-2)/24)-b4*P4*((dt(j))^(-lamda4))*((((-h(j+1)*dx)/q(j+1))-2)/24);
            
            A(i,i)   = ((2*q(j+1))/(dx)^2)-((1/6)*(g(i,j+1)-(((h(j+1))^2)/(q(j+1))^2)))+(g(i,j+1)+(((dx)^2)/12)*(g2(i,j+1)-((h(j+1)*g1(i,j+1))/(q(j+1)))))-b1*P1*((dt(j))^(-lamda1))*(-5/6)-b2*P2*((dt(j))^(-lamda2))*(-5/6)-b3*P3*((dt(j))^(-lamda3))*(-5/6)-b4*P4*((dt(j))^(-lamda4))*(-5/6);
            
            A(i,i+1) = (-q(j+1)/(dx)^2)+(1/12)*(g(i,j+1)-(((h(j+1))^2)/(q(j+1))^2))+(1/2*dx)*(h(j+1)+((dx)/24)*(2*(g1(i,j+1))-((h(j+1)*g(i,j+1))/q(j+1))))-b1*P1*((dt(j))^(-lamda1))*((((h(j+1)*dx)/q(j+1))-2)/24)-b2*P2*((dt(j))^(-lamda2))*((((h(j+1)*dx)/q(j+1))-2)/24)-b3*P3*((dt(j))^(-lamda3))*((((h(j+1)*dx)/q(j+1))-2)/24)-b4*P4*((dt(j))^(-lamda4))*((((h(j+1)*dx)/q(j+1))-2)/24);
    end
    
    %m+1th row elements by using the u(1,t)=1+t^(alpha)
    A(M+1,M+1)=1;
    
    for i=2:M
        for t1=1:j-1
            sum1(i,j)=sum1(i,j)+((dt(t1))^(-1))*P1*((t(j+1)-t(t1))^(1-lamda1)-(t(j+1)-t(t1+1))^(1-lamda1))*(((((-h(j+1)*dx)/q(j+1))-2)/24)*(u(i-1,t1+1)-u(i-1,t1))...
                +(-5/6)*(u(i,t1+1)-u(i,t1))+((((h(j+1)*dx)/q(j+1))-2)/24)*(u(i+1,t1+1)-u(i+1,t1)));
        end
    end
   
     for i=2:M
        for t1=1:j-1
            sum2(i,j)=sum2(i,j)+((dt(t1))^(-1))*P2*((t(j+1)-t(t1))^(1-lamda2)-(t(j+1)-t(t1+1))^(1-lamda2))*(((((-h(j+1)*dx)/q(j+1))-2)/24)*(u(i-1,t1+1)-u(i-1,t1))...
                +(-5/6)*(u(i,t1+1)-u(i,t1))+((((h(j+1)*dx)/q(j+1))-2)/24)*(u(i+1,t1+1)-u(i+1,t1)));
        end
     end

     for i=2:M
        for t1=1:j-1
            sum3(i,j)=sum3(i,j)+((dt(t1))^(-1))*P3*((t(j+1)-t(t1))^(1-lamda3)-(t(j+1)-t(t1+1))^(1-lamda3))*(((((-h(j+1)*dx)/q(j+1))-2)/24)*(u(i-1,t1+1)-u(i-1,t1))...
                +(-5/6)*(u(i,t1+1)-u(i,t1))+((((h(j+1)*dx)/q(j+1))-2)/24)*(u(i+1,t1+1)-u(i+1,t1)));
        end
     end

     for i=2:M
        for t1=1:j-1
            sum4(i,j)=sum4(i,j)+((dt(t1))^(-1))*P4*((t(j+1)-t(t1))^(1-lamda4)-(t(j+1)-t(t1+1))^(1-lamda4))*(((((-h(j+1)*dx)/q(j+1))-2)/24)*(u(i-1,t1+1)-u(i-1,t1))...
                +(-5/6)*(u(i,t1+1)-u(i,t1))+((((h(j+1)*dx)/q(j+1))-2)/24)*(u(i+1,t1+1)-u(i+1,t1)));
        end
     end

    B(1,1)=u(1,j);
    for i=2:M

           B(i,1)=b1*sum1(i,j)+b2*sum2(i,j)+b3*sum3(i,j)+b4*sum4(i,j)-b1*P1*((dt(j))^(-lamda1))*(((((-h(j+1)*dx)/q(j+1))-2)/24)*u(i-1,j)+(-5/6)*u(i,j)+((((h(j+1)*dx)/q(j+1))-2)/24)*u(i+1,j))...
                 -b2*P2*((dt(j))^(-lamda2))*(((((-h(j+1)*dx)/q(j+1))-2)/24)*u(i-1,j)+(-5/6)*u(i,j)+((((h(j+1)*dx)/q(j+1))-2)/24)*u(i+1,j))...
                 -b3*P3*((dt(j))^(-lamda3))*(((((-h(j+1)*dx)/q(j+1))-2)/24)*u(i-1,j)+(-5/6)*u(i,j)+((((h(j+1)*dx)/q(j+1))-2)/24)*u(i+1,j))...
                 -b4*P4*((dt(j))^(-lamda4))*(((((-h(j+1)*dx)/q(j+1))-2)/24)*u(i-1,j)+(-5/6)*u(i,j)+((((h(j+1)*dx)/q(j+1))-2)/24)*u(i+1,j))...
                 +((((h(j+1)*dx)/q(j+1))+2)/24)*s(i-1,j+1)+(5/6)*s(i,j+1)+((((-h(j+1)*dx)/q(j+1))+2)/24)*s(i+1,j+1);
  
    end
    
    B(M+1,1)=u(M+1,j);
    
    r1=A\B;
   
    for i=2:M
        u(i,j+1)=r1(i);
    end 
end

 for i=1:M+1
    for j=1:N+1
        err1(i,j)=(e(i,j)-u(i,j))^2;
    end
end


 for i=1:M+1
   for j=1:N+1
       err(i,j)=abs(e(i,j)-u(i,j));
   end
 end

 %E1=e(6,:);
 %E2=u(6,:);
err1=err(:,N+1);
err2=err1';
err3=(err1).^2;
sum4=0;
for i=1:M+1
     sum4=sum4+err3(i);                    % L2-err in space
end 
L2=sqrt(dx*sum4)

% mae=max(err(:))                         % L_inf error

%sum5=0;
%for i=2:M
%    sum5=sum5+err1(i);                            % L1-err
%end
%L1=dx*sum5


%err2 = sum(err1);
% disp(err1)
% disp(err2)
%err3 = max(err2);
%L2 = sqrt(dx*err3)                %L2 error in time

%plot(w,t,'-^','Color','r','Linewidth',2)
% subplot(3,1,1)
% surf(x,t,transpose(err1))
% subplot(3,2,1)
% surf(x,t,e)

%save('surface_error_graded_04.mat','err')
%save('surface_error_uniform_08.mat','err')
%save('surface_exact_uniform_08.mat','e') 
%save('surf_computed_graded_08.mat','u') 
%save time04.m t -ascii -double
%save x.m x -ascii -double
%save t_graded_08.m t -ascii -double
%save x_graded.m x -ascii -double

% Ts=t';
% t1=t(2:2:end) %even
% x1=x(2:2:end) %even
% t2=t(1:2:end) %odd
% save timegrid.m w -ascii -double
 %save exact_uniform_04.m E1 -ascii -double
%save time_uniform.m t -ascii -double
% save timegrid04.m t -ascii -double
%  save uniform_exact_03.m E1 -ascii -double
%  save approximate_graded_03.m E2 -ascii -double
% save timetime.m t -ascii -double
toc