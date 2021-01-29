% Propagador de Fourier
clear all; close all;
set(0,'defaultTextInterpreter','latex');

N = 2^10;                % Numero de puntos, escogemos una potencia de 2
                        % para que el calculo de la transformada sea más
                        % rápido

% Creamos el vector de indices
NV = (-N/2:1:N/2-1);     % Con esto nos aseguramos de tener la frecuencia 0
                         % y de que NV tiene N elementos
L = 1e-3;                % Las unidades son metros
dx = 2*L/N;
kmax = pi/dx;
[X,Y] = meshgrid(NV*dx);
[Kx,Ky] = meshgrid(kmax*2/N*NV);

lambda = 633e-9;
k = 2*pi/lambda;

r = sqrt(X.^2+Y.^2);
kt = sqrt(Kx.^2+Ky.^2);

w0 = 0.15e-3;                            % Cintura del haz Gaussiano
zR = k*w0^2/2;
z = 1.5*zR;      %Distancia de propagación maxima en unidades de la distancia de Rayleigh
nz = 300;
dz = z/nz;

%****************************** Propagador ******************************
Prop = (exp(-1i*0.5*dz*(kt.^2)/k));      %Full paraxial propagator
%Prop = exp(i*dz*sqrt(k^2-kt.^2));       %Non paraxial propagator

%************************* Perfil Inicial *******************************
flens = 0.05; %0.3*zR;
Tlens = 1;
Tlens = exp(-1i*k/(2*flens)*r.^2);             % Funcion de transmitancia de lente
f = exp(-r.^2/w0^2).*Tlens;                    % Haz Gaussiano

Ur = zeros(N,nz+1);
U0 = f;
Ur(:,1) = U0(:,N/2+1);

%******************* Inicia la propagacion ******************************
F = fftshift(fft2(U0));
for ii=1:nz
    F = F.*Prop;
    A = ifft2(F);
    Ur(:,ii+1)=A(:,N/2+1);
end
%Ur = ifft2(Ur);

%*******************  Grafica del campo inicial *************************
figure(2),subplot(1,2,1),surf(X(N/2+1,:)/w0,Y(:,N/2+1)/w0,abs(A));shading interp,lighting phong, view(2)
colormap parula;
ejes = gca;
ejes.FontSize = 13;
title('Campo propagado con Fourier')
xlabel('$x/w_{0}$','FontSize',20);
ylabel('$y/w_{0}$','FontSize',20);
axis square;


%*******************  Grafica del campo analítico ***********************
U = exp(1i*k*L)/(1i*(1.5-1i*(1-30*zR))).*exp(r.^2./zR.*k./3.*(1i-1/(1.5-1i*(1-30*zR))));
figure(2),subplot(1,2,2),surf(X(N/2+1,:)/w0,Y(:,N/2+1)/w0,abs(U));shading interp,lighting phong, view(2)
colormap parula; ejes = gca;
ejes.FontSize = 13;
title('Campo Analitico')
xlabel('$x/w_{0}$','FontSize',20);
ylabel('$y/w_{0}$','FontSize',20);
axis square;
%%
figure(2);
surf((0:dz:z)/zR,Y(:,N/2+1)/w0,abs(Ur)),shading interp,colormap parula, view(2);   
ejes = gca; axis([0 1.5 -6 6])
ejes.FontSize = 13;
title('Propagacion del campo con lente')
xlabel('$z/z_{R}$','FontSize',20);
ylabel('$y/w_{0}$','FontSize',20);
axis square