function kf = testeDM(A,B,C,K,Q,R,Pf,T,N,x0,r,m,M,g,L,Sx,bx,Su,bu,Sf,bf,Oinf)
%% Função teste para o controlador dual mode
% Os parâmetros de entrada são:
% kf = testeDM(A,B,C,K,Q,R,Pf,T,N,x0,r,m,M,g,L,Sx,bx,Su,bu,Sf,bf,Oinf)
% Do modelo xp = Ax + Bu discretizado:
% A  ----------- Matriz do modelo discretizado A
% B  ----------- Matriz do modelo discretizado B
% C  ----------- Matriz do modelo discretizado C
% K  ----------- Matriz do regulador Dual Mode
% Q  ----------- Matriz de pesos do estado
% R  ----------- Matriz de pesos do controle
% Pf ----------- Matriz solução da equação de Ricatti, DLQR
% T  ----------- Período de amostragem adotado
% N  ----------- Número de iterações para o horizonte móvel
% x0 ----------- Estado inicial
% r  ----------- Referência (estado final desejado)
% m  ----------- Massa do pêndulo
% M  ----------- Massa do carrinho
% g  ----------- Aceleração da gravidade
% L  ----------- Comprimento do pêndulo
% Sx ----------- Matriz de limites do estado Sx x <= bx
% bx ----------- Matriz de limites do estado Sx x <= bx
% Su ----------- Matriz de limites do controle Su u <= bu
% bu ----------- Matriz de limites do controle Su u <= bu
% Sf ----------- Matriz de limites Oinf
% bf ----------- Matriz de limites Oinf
% Oinf --------- Poliedro de limites, no R^4

n = size(A,1);
p = 1;
q = size(C,1);
aux = [A-eye(n) B; C zeros(q,1)]\[zeros(n,q);eye(q)];
Nx = aux(1:n,:);
Nu = aux(n+1:end,:);
xbar = Nx*r;
ubar = Nu*r;

equacao_estado = @(t,x,u) [x(2); m/M*x(3) - L*x(4)^2*x(3) + u/M; x(4);...
                           (m/M+1)*g/L*x(3) - x(4)^2*x(3) + u/(M*L)];                      

x{1} = x0;                       
xs{1} = x0;
trec = 0; xrec = xs{1}';
k = 1;
while(max(Sf*x{k} - bf) > 0)
    % Controle
    un = otimizador_u(A,B,Q,R,Pf,N,Sx,bx,Su,bu,Sf,bf,x{k},xbar,ubar);
    u(k) = un(1);
    x{k+1} = A*x{k} + B*u(k);
    % Simulação a tempo "contínuo"
    delta_x{k} = xs{k} - x0;
    uns = otimizador_u(A,B,Q,R,Pf,N,Sx,bx,Su,bu,Sf,bf,delta_x{k},xbar,ubar);
    delta_u(k) = uns(1);
    us(k) = delta_u(k);
    [tout, xout] = ode45(@(t,x) equacao_estado(t,x,us(k)),[(k-1)*T k*T], xs{k});
    trec = [trec; tout(2:end)];
    xrec = [xrec; xout(2:end,:)];
    xs{k+1} = xrec(end,:)';
    k = k+1;
end
kf = k;
for k = kf:kf+30
    % Controle
    u(k) = ubar - K*(x{k} - xbar);
    x{k+1} = A*x{k} + B*u(k);
    % Simulação a tempo "contínuo"
    delta_x{k} = xs{k} - x0;
    delta_u(k) = -K*(delta_x{k} - xbar) + ubar;
    us(k) = delta_u(k);
    [tout, xout] = ode45(@(t,x) equacao_estado(t,x,us(k)),[(k-1)*T k*T], xs{k});
    trec = [trec; tout(2:end)];
    xrec = [xrec; xout(2:end,:)];
    xs{k+1} = xrec(end,:)';
end
kf = k;

DT1 = delaunayTriangulation(Oinf.V(:,1:2)); 
DT2 = delaunayTriangulation(Oinf.V(:,3:4));
index1 = convexHull(DT1);
index2 = convexHull(DT2);
y = cell2mat(x);

figure;
hold on;
p1 = patch(DT1.Points(index1,1),DT1.Points(index1,2),'r','DisplayName',"$\chi_{f,\perp d\ vs\ \dot{d}}$",'LineStyle','none');
p2 = plot(xrec(:,1),xrec(:,2),'k-','linewidth',2,"DisplayName","Cont\'{i}nuo");
p3 = plot(y(1,:),y(2,:),'b-o','linewidth',1,'DisplayName',"Dual Mode");
p4 = plot(r(1),0,'yo','markersize',5,'linewidth',5,"DisplayName","Objetivo");
ylabel('$\dot{d}\ [m/s]$','Interpreter','latex', 'FontSize', 12);
xlabel('$d\ [m]$', 'FontSize', 12,'Interpreter','latex');
%title("Gr\'{a}fico $d(t)$ x $\dot{d}(t)$", 'FontSize', 14,'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
grid on;
hleg = legend([p2 p3 p1 p4],'location','best','FontSize', 12,'Interpreter','latex');
print(gcf,'/Users/Reynaldo/Documents/ITA/EE-254/figures/dm_1.eps','-depsc','-r300');  
print(gcf,'/Users/Reynaldo/Documents/ITA/EE-254/figures/dm_1.png','-dpng','-r300');  

figure;
hold on;
p1 = patch(DT2.Points(index2,1),DT2.Points(index2,2),'r','DisplayName',"$\chi_{f,\perp \theta\ vs\ \dot{\theta}}$",'LineStyle','none');
p2 = plot(xrec(:,3),xrec(:,4),'k-','linewidth',2,"DisplayName","Cont\'{i}nuo");
p3 = plot(y(3,:),y(4,:),'b-o','linewidth',1,"DisplayName","Dual Mode");
p4 = plot(r(2),0,'yo','markersize',5,'linewidth',5,'DisplayName',"Objetivo");
ylabel('$\dot{\theta}\ [rad/s]$','Interpreter','latex', 'FontSize', 12);
xlabel('$\theta\ [rad]$', 'FontSize', 12,'Interpreter','latex');
%title("Gr\'{a}fico $\theta(t)$ x $\dot{\theta}(t)$", 'FontSize', 14,'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
grid on;
hleg = legend([p2 p3 p1 p4],'location','best','FontSize', 12,'Interpreter','latex');
print(gcf,'/Users/Reynaldo/Documents/ITA/EE-254/figures/dm_2.eps','-depsc','-r300');  
print(gcf,'/Users/Reynaldo/Documents/ITA/EE-254/figures/dm_2.png','-dpng','-r300');  

hold off;
figure;
p1 = stairs(0:T:(k-1)*T,us,'k-','linewidth',2,'DisplayName',"Cont\'{i}nuo");
hold on;
p2 = stairs(0:T:(k-1)*T,u,'b-','linewidth',1,'DisplayName',"Dual Mode");
p3 = yline(-bu(2)/Su(1),'r--','linewidth',2,'DisplayName',"Limites");
p4 = yline(-bu(1)/Su(2),'r--','linewidth',2,'DisplayName',"Limites");
p5 = yline(0,'y','linewidth',2,'DisplayName',"Objetivo");
xlabel('$t\ [s]$','Interpreter','latex', 'FontSize', 12);
ylabel('$u\ [N]$', 'FontSize', 12,'Interpreter','latex');
%title("Gr\'{a}fico $u(t)$ x $t$", 'FontSize', 14,'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',[0 (kf)*T]);
grid on;
hleg = legend([p1 p2 p3 p5],'location','best','FontSize', 12,'Interpreter','latex');
print(gcf,'/Users/Reynaldo/Documents/ITA/EE-254/figures/dm_3.eps','-depsc','-r300');  
print(gcf,'/Users/Reynaldo/Documents/ITA/EE-254/figures/dm_3.png','-dpng','-r300');  


hold off;
figure;
p1 = plot(trec,xrec(:,3),'k-','linewidth',2,'DisplayName',"Cont\'{i}nuo");
hold on;
p2 = plot(0:T:(kf)*T,y(3,:),'b-o','linewidth',1,'DisplayName',"Dual Mode");
p3 = yline(-bx(2)/Sx(1,3),'r--','linewidth',2,'DisplayName',"Limites");
p4 = yline(-bx(1)/Sx(2,3),'r--','linewidth',2,'DisplayName',"Limites");
p5 = yline(0,'y','linewidth',2,'DisplayName',"Objetivo");
xlabel('$t\ [s]$','Interpreter','latex', 'FontSize', 12);
ylabel('$\theta\ [rad]$', 'FontSize', 12,'Interpreter','latex');
%title("Gr\'{a}fico $\theta(t)$ x $t$", 'FontSize', 14,'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','XLim',[0 (kf)*T]);
grid on;
hleg = legend([p1 p2 p3 p5],'location','best','FontSize', 12,'Interpreter','latex');
print(gcf,'/Users/Reynaldo/Documents/ITA/EE-254/figures/dm_4.eps','-depsc','-r300');  
print(gcf,'/Users/Reynaldo/Documents/ITA/EE-254/figures/dm_4.png','-dpng','-r300');  

close all;
end