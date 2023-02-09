clear; clc; close all;
%% Problem 2(c)
a = 4; b = 2;
delta = pi/3;
n = 100;

t = linspace(0,2*pi,n);
x = sin(a*t*delta);
y = sin(b*t);

u = zeros(1,n-2);
v = zeros(1,n-2);
xy = [x; y];
diff = xy(:,2:end)-xy(:,1:end-1);

%%% YOUR CODE TO COMPUTE GRADIENT HERE %%%
%%% END HOMEWORK PROBLEM %%%

figure;
plot(x,y,'linewidth',2,'color','black'); hold on;
quiver(x(2:end-1),y(2:end-1),u,v,'linewidth',1,'color','red');
axis equal;

%% Problem 2(d)
%%% YOUR CODE TO COMPUTE KAPPA HERE %%%
%%% END HOMEWORK PROBLEM %%%

figure;
X = x(2:end-1);
Y = y(2:end-1);
Z = zeros(size(kappa));
surface([X;X],[Y;Y],[Z;Z],[kappa;kappa],...
    'edgecolor', 'interp', 'linewidth',2);
axis equal;
colorbar;

%% Problem 2(e)
t0 = 0;
t1 = pi*1.25;
nSamples = 100;
nSteps = 20;

% We provide a few examples of curves to try
%curveFunction = @(t) [(cos(t)-cos(3*t).^3); (sin(t)-sin(3*t).^3)]';
%curveFunction = @(t) [cos(t);sin(t)]';
curveFunction = @(t) [t;(t-t0).*(t1-t)]';
curve = curveFunction(linspace(t0,t1,nSamples));

% Time step
f = figure;
plt = plot(curve(:,1),curve(:,2),'k','linewidth',2);
axis equal;
for i=1:nSteps
    %%% YOUR CODE HERE TO PERFORM GRADIENT DESCENT %%%
    %%% END HOMEWORK PROBLEM %%%
    plt.XData = curve(:,1);
    plt.YData = curve(:,2);
    drawnow;
end

