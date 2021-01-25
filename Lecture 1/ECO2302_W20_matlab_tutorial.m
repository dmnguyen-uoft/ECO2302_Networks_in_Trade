%% commenting
% This file provides a basic introduction to some of the key Matlab
% commands.
%{
  Comment 1.
  Comment 2.
  Comment 3.
%}

%% clearing workspace, figures, command line
clearvars
clear global
close all
clc

%% getting help
%help sqrt

%% matrices, cell arrays, and data structures

% entering matrices manually
A = [1 2 3;
     4 5 6; 
     7 8 9];

% building number sequences
B = 0:2:10;
C = linspace(1,100,5);

% zeros and ones matrices
D = zeros(3,4);
E = ones(3,4);

% concatenating matrices
F = [D;E];

% accessing matrix elements
a = A(3,2);
arow = A(1,:);
acol = A(:,1);
b = B(end);

% dimensions
A_size = size(A);
D_size = size(D);

% cell arrays
G = cell(3,1);
G{1} = 10;
G{2} = 'abc';
G{3} = @(x) 2*exp(x);

% structures
H.first_name = 'Justin';
H.last_name = 'Trudeau';
H.age = 47;

%% matrix operations and functions

% addition and subtraction
Aplus = A + ones(size(A));
Aminus = A - ones(size(A));

% multiplication
Atimes = A * ones(size(A'));    % matrix product
Adtimes = A .* ones(size(A));   % element by element product

% division
K1 = [1 .1 .2;
     .1 1 .3;
     .2 .3 1];
K2 = [1 .2 .3;
      .5 1 .1;
      .4 .1 1];
Kdivr = K1/K2;  % same as K1 * inv(K2), but faster
Kdivl = K1\K2;  % same as inv(K1) * K2, but faster
Kdivd = K1./K2; % element by element division

% exponentiating
Apower = A^2;   % same as A * A
Adpower = A.^2; % element by element exponent

% inverse
K = [1 .1 .2;
     .1 1 .3;
     .2 .3 1];
Kinv = K^-1;    % same as inv(K)

% column and row sums
Kcolsum = sum(K);   % same as sum(K,1)
Krowsum = sum(K,2);

% maximum and minimum elements
[Kmax,Imax] = max(K);    % Kmax is maximum value in each column of K; Imax is row index of maximum value in each column of K
[Kmin,Imin] = min(K);    % Kmin is minimum value in each column of K; Imin is row index of minimum value in each column of K

% sorting
[Ksort,Isort] = sort(K(:)); % Ksort is sorted values of K in ascending order; Isort is sorted index of each element of K

% eigenvalue decomposition
%   eigVec is a matrix where each column is an eigenvector of K
%   eigVal is diagonal matrix where each diagonal element is an eigen value of K
[eigVec,eigVal] = eig(K);

%% random variables

x = rand;           % uniform RV on [0,1] same as unifrnd(0,1)

x = normrnd(2,1);   % normal RV with mean 2 and s.d. 1

x = normpdf(0,2,1); % pdf of normal RV with mean 2 and s.d. 1, evaluated at 0
x = normcdf(0,2,1); % cdf of normal RV with mean 2 and s.d. 1, evaluated at 0


%% functions

% anonymous functions
f = @(x,y,z) x^2 + .5 * y^3 - sqrt(z);
fval1 = f(1.3,2.9,8);

% function scripts (see file funscript.m)
out = funscript(1.3,2.9,8);

%% if/for/while loops

% if statements
%disp(x)
if x > 0
    y = 100;
elseif x < 0
    y = 200;
else
    y = 300;
end

% for loops
N = 4;
zvec = zeros(N,1);
nvec = [1 4 6 10];

for i = 1:length(nvec)
    ni = nvec(i);
    zvec(i) = ni^2;
    
end

% while loops
eps = 1;
while eps>.01
    eps = eps*.9; 
end


%% plotting: example 1 - sine and cosine waves

% define vector of x-axis values for plotting
xvec = linspace(-2*pi,2*pi,100);

% initialize figure
figure();

% set 'hold' to on (this lets you plot multiple things in the same figure)
hold on     

% plot sine wave
plot(xvec,sin(xvec),'linewidth',2)  

% plot cosine wave
plot(xvec,cos(xvec),'linewidth',2,'linestyle',':')    

% set figure properties
xlabel('x')                         % x-axis label
ylabel('value')                     % y-axis label
title('sine and cosine waves')      % title
legend('sin(x)','cos(x)')           % legend
set(gca,'fontsize',20)              % font size


%% plotting: example 2 - Batman

% define functions for plotting
f1 = @(x) 1.5*sqrt((-abs(abs(x)-1)).*abs(3-abs(x))./((abs(x)-1).*(3-abs(x)))).*(1+abs(abs(x)-3)./(abs(x)-3)).*sqrt(1-(x./7).^2)+(4.5+.75.*(abs(x-.5)+abs(x+.5))-2.75.*(abs(x-.75)+abs(x+.75))).*(1+abs(1-abs(x))./(1-abs(x)));
f2 = @(x) (-3).*sqrt(1-(x./7).^2).*sqrt(abs(abs(x)-4)./(abs(x)-4));
f3 = @(x) abs(x./2) - .0913722.*x.^2 -3 + sqrt(1-(abs(abs(x)-2)-1).^2);
f4 = @(x) (2.71052+1.5-.5.*abs(x)-1.35526.*sqrt(4-(abs(x)-1).^2)).*sqrt(abs(abs(x)-1)./(abs(x)-1));

% define x-axis values for plotting
xvec1 = [-7:.01:-3 -1:.0001:1 3:.01:7];
xvec2 = [-7:.01:-4 4:.01:7];
xvec3 = -4:.01:4;
xvec4 = [-3:.01:-1 1:.01:3];

% initialize figure and figure property values
figure(); lwidth = 2; fsize = 26;

% set 'hold' to on (this lets you plot multiple things in the same figure)
hold on  

% plot each part of the figure
plot(xvec1,f1(xvec1),'LineWidth',lwidth,'Color','k','LineStyle','-')
plot(xvec2,f2(xvec2),'LineWidth',lwidth,'Color','b','LineStyle','--')
plot(xvec3,f3(xvec3),'LineWidth',lwidth,'Color','r','LineStyle',':')
plot(xvec4,f4(xvec4),'LineWidth',lwidth,'Color','g','LineStyle','-.')

% set figure properties
xlabel('x','FontSize',fsize)    % x-axis label
ylabel('y','FontSize',fsize)    % y-axis label
title('Batman')                 % title
legend('upper wings and head','lower wings','lower body','neck')    % legend
set(gca,'FontSize',fsize,'Xlim',[-8 13])    % font size and x-axis limits

%% numerical solvers

% set up equation to be solved
xvec = linspace(.01,1,100);
f = @(x) log(x) + exp(x) - 1;
figure(); hold on
plot(xvec,f(xvec))
plot([0 1],[0 0])

% fsolve
options = optimoptions('fsolve','display','iter');
xinit = .1;
xsol = fsolve(f,xinit,options);

%% numerical optimization

% set up function to be minimized
f = @(x) (x(1)-2)^2 + (x(2)-3)^2;
%f = @(x) (x(1)-2)^2 + (x(2)-3)^2 + log(x(1) + x(2));
Nx = 100;
x1 = linspace(.01,5,Nx);
x2 = linspace(.01,5,Nx);
fvals = zeros(Nx);
for i = 1:Nx
   for j = 1:Nx
       fvals(j,i) = f([x1(i),x2(j)]);
   end
end
figure(); 
subplot(121); surf(x1,x2,fvals)
subplot(122); contour(x1,x2,fvals)

% fminsearch
options = optimset('display','iter');
xinit = [1;1];
xmin = fminsearch(f,xinit,options);

