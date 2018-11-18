close all
clear all

star_a = 0.8;
star_b = 0.3;
star_k = 5.0;

phi = 0:0.01:2*pi;
r   = star_a + star_b*sin(star_k*phi);
x   = r.*cos(phi);
y   = r.*sin(phi);

plot(x,y)
legend({'$r = 0.5 + 0.6*\sin(2.0*\phi), 0 \leq \phi \leq 2\pi$'},'Interpreter','latex')
xlabel('x-axis','fontsize',20);
ylabel('y-axis','fontsize',20);
axis([-1.99 1.99 -1.99 1.99]);
axis('equal');

hold on

