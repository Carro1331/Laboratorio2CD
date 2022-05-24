clc;
clear;
close all;
Fm = 1000 ;
Tm = 1/Fm;
t = (-0.05:Tm:0.05);
t2 = (-0.5:Tm:0.5);
%%%%%%%%% GRAFICO RESPUESTA AL IMPULSO COSENO ALZADO DE NYQUIST %%%%%%%%%%%
figure(1)
%subplot(2,1,1);
hold on
grid on
plot(rollof(0));plot(rollof(0.25));plot(rollof(0.75));plot(rollof(1));
xlim([400 600]);ylim([-5 22]);

title('RESPUESTA AL IMPULSO COSENO ALZADO DE NYQUIST'); xlabel('Tiempo(s)'); ylabel('AMPLITUD(m)');
legend('rolloff: 0','rolloff: 0.25','rolloff: 0.75','rolloff: 1');
hold off

%%%%%%% GRAFICO RESPUESTA EN FRECUENCIA COSENO ALZADO DE NYQUIST %%%%%%%%%%
figure(2)
%subplot(2,1,1);
hold on
[f,FFT] = TransformadaFourier(rollof(0),1000,Fm);
[f1,FFT1] = TransformadaFourier(rollof(0.25),1000,Fm);
[f2,FFT2] = TransformadaFourier(rollof(0.75),1000,Fm);
[f3,FFT3] = TransformadaFourier(rollof(1),1000,Fm);
plot(f, FFT);plot(f1, FFT1);plot(f2, FFT2);plot(f3, FFT3);
xlim([-200 200]);ylim([0 1.3]);
xlabel('FRECUENCIA(Hz)');ylabel('AMPLITUD(m)');
title('RESPUESTA EN FRECUENCIA COSENO ALZADO DE NYQUIST');
legend('roll-of:0','roll-of:0.25','roll-of:0.75','roll-of:1');

%%%%%%%%%%%%%% DIAGRAMA DE OJO COSENO ALZADO DE NYQUIST %%%%%%%%%%%%%%%%%%%
figure(3)
newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};
colororder(newcolors);

for j=1:1000 %%SE PUEDE REDUCIR ESTE VALOR PARA MEJORAR LA EJECUCION DE CODIGO 

    s = 2*randi([0,1],1,3)-1; %Generamos los numeros aleatorios
    s = awgn(s,20,'measured'); %Agregamos el ruido gauseano
    pulseT = upsample(s, length(t)); %Agrega muestras
    
    aux = pulseT(1);
    sizeX = length(pulseT);
    sizet =length(t);
    c = 1;
    for i=1:sizeX
        if c*sizet+1 == i
            c = c+1;
            aux = pulseT(i);
        end
        pulseT(i) = aux;
    end  
    
    finalPulse = conv(pulseT,rollof2(0.22)); %aplicamos, la convolucion debida
    hold on
    subplot(2,1,1);
    plot(finalPulse);
    %xlim([100 300]);ylim([-1000 1000 ]); %Se aplica el corte necesario para mostrar el diagrama de ojo en la mejor manera.
end

title('DIAGRAMA DE OJO DE COSENO ALZADO DE NYQUIST CON ROLLOFF DE 0.22');
xlabel('TIEMPO(s)');ylabel('AMPLITUD(m)');
hold off
grid on
subplot(2,1,2);
plot(rollof2(0.22));
xlim([0 100]);ylim([-5 22]);
xlabel('TIEMPO(s)');ylabel('AMPLITUD(m)');
legend('rolloff:0.22');
title('COSENO ALZADO DE NYQUIST CON ROLLOFF DE 0.22');


%%%%%% Coseno Alzado de Nyquis, aplicado a un Rolloff solicitado. %%%%%%%%

%%%%%Existen dos funciones exactamente iguales, pero trabajan con distintos
%%%%%intervalos de tiempo, para permitir una mejor visualizacion y trabajo
%%%%%de la señal.
function cosalzado = rollof(r) 
    Fm = 1000;          % Frecuencia de Muestreo
    Tm = 1/Fm;          % Tiempo de Muestreo
    t = (-0.5:Tm:0.5);  % Vector de Tiempo
    fo = 10;            % Ancho de Banda del filtro
    B = fo*(r+1);       % Ancho de Banda absoluto 
    fd = B-fo;          % Parámetro fd

    % PARTE UNO
    numsen = 2*fo*sin(2*fo*pi*t);          % Numerador de la función impulso de coseno alzado (Seno).
    densen = (2*fo*pi*t);                  % Denominador de la función impulso de coseno alzado (Seno).
    cerossen = find(abs(densen) == 0);     % Encontrar los ceros del denominador.
    operarsen = numsen./densen;            % Primera operación de la función impulso de coseno alzado (Seno).
    operarsen(cerossen) = 2*fo;            % Donde se encontraron los ceros del denominador, se iguala al valor de 2*fo.

    % PARTE DOS
    numcos = cos(2*fd*pi*t);               % Numerador de la función impulso de coseno alzado (Coseno).
    dencos = (1-(4*fd*t).^2);              % Denominador de la función impulso de coseno alzado (Coseno).
    ceroscos = find(abs(dencos) == 0);     % Encontrar los ceros del denominador.
    operarcos = numcos./dencos;            % Segunda operación de la función impulso de coseno alzado (Coseno).
    operarcos(ceroscos) = pi/4;            % Donde se encontraron los ceros del denominador, se iguala al valor de pi/4.

    cosalzado = operarsen.*operarcos;      % Operación final para tener la función impulso de coseno alzado de nyquist.
end

function cosalzado = rollof2(r) 
    Fm = 1000;          % Frecuencia de Muestreo
    Tm = 1/Fm;          % Tiempo de Muestreo
    t = (-0.05:Tm:0.05);  % Vector de Tiempo
    fo = 10;            % Ancho de Banda del filtro
    B = fo*(r+1);       % Ancho de Banda absoluto 
    fd = B-fo;          % Parámetro fd

    % PARTE UNO
    numsen = 2*fo*sin(2*fo*pi*t);          % Numerador de la función impulso de coseno alzado (Seno).
    densen = (2*fo*pi*t);                  % Denominador de la función impulso de coseno alzado (Seno).
    cerossen = find(abs(densen) == 0);     % Encontrar los ceros del denominador.
    operarsen = numsen./densen;            % Primera operación de la función impulso de coseno alzado (Seno).
    operarsen(cerossen) = 2*fo;            % Donde se encontraron los ceros del denominador, se iguala al valor de 2*fo.

    % PARTE DOS
    numcos = cos(2*fd*pi*t);               % Numerador de la función impulso de coseno alzado (Coseno).
    dencos = (1-(4*fd*t).^2);              % Denominador de la función impulso de coseno alzado (Coseno).
    ceroscos = find(abs(dencos) == 0);     % Encontrar los ceros del denominador.
    operarcos = numcos./dencos;            % Segunda operación de la función impulso de coseno alzado (Coseno).
    operarcos(ceroscos) = pi/4;            % Donde se encontraron los ceros del denominador, se iguala al valor de pi/4.

    cosalzado = operarsen.*operarcos;      % Operación final para tener la función impulso de coseno alzado de nyquist.
end

%%% Transformada de Fourier para el impulso de coseno alzado de nyquist %%%
function [f, FFT] = TransformadaFourier(s, fs, Fm)
    FFT = abs(fftshift(fft(s, 1001)))/fs;
    f = [-500:500]*(Fm/100);
end
