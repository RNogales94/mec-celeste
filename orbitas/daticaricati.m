%Datos sobre los Planetas
%Planeta = [Duración del año en días terrestres, i, OMEGA, a, epsilon, omega*, T]

clear; clc;

mi=398600;
scarto=1.e-8;

Mercurio=[87.97, 7.00, 47.14, 0.387, 0.206, 75.90];
Venus=[224.7, 3.59, 73.78, 0.723, 0.007, 130.15];
Tierra=[365.26, 0.00, 0.00, 1.000, 0.017, 101.22];
Marte=[686.98, 1.85, 48.78, 1.524, 0.093, 334.22];
Jupiter=[4332.6, 1.31, 99.44, 5.203, 0.048, 12.72];
Saturno=[10759, 2.5, 112.79, 9.546, 0.056, 91.09];
Urano=[30687, 0.77, 73.48, 19.20, 0.047, 169.05];
Neptuno=[60784, 1.78, 130.68, 30.09, 0.009, 43.83];

%Elegir un planeta y un instante de tiempo. La persona que quiere usar el
%programa tiene que elegir un Planeta y el instante de tiempo, el programa
%dice que posicion tiene el planeta elegido en el instante de tiempo
%elegido.
%P es el planeta, t el instante de tiempo, m es necesario para graficar la
%elipse .

P=input('Elegir el Planeta: ');

t=input('Elegir un instante de tiempo: ');

m=1000;

l=linspace(0,P(1),m);

%Usamos el metodo de Newton para encontrar la posicion del planeta en el
%instante t (y en otros instantes, para graficar la elipse).
i=1;
while i<=m
    
    T=(2*pi*P(4)^(1.5)/sqrt(mi))/24*60*60;

    M=2*pi*l(i)/T;

    if M<pi
        u0=M+P(5)/2;
    else
        u0=M-P(5)/2;
    end

    r=1;

    while abs(r)>scarto
        r=(u0-P(5)*sin(u0)-M)/(1-P(5)*cos(u0));
        u0=u0-r;
    end

    fx(i)=P(4)*[cos(u0)-P(5)];
    sx(i)=P(4)*[sqrt(1-P(5)^2)*sin(u0)];
    i=i+1;
end

%Periodo orbital, T
T=(2*pi*P(4)^(1.5)/sqrt(mi))/24*60*60;

%Anomalia media
M=2*pi*t/T;

if M<pi
    u0=M+P(5)/2;
else
    u0=M-P(5)/2;
end

r=1;
while abs(r)>scarto
    r=(u0-P(5)*sin(u0)-M)/(1-P(5)*cos(u0));
    u0=u0-r;
end

u0
x=P(4)*[cos(u0)-P(5), sqrt(1-P(5)^2)*sin(u0)]

plot(fx,sx,'r*');
hold on
plot(x(1),x(2),'b*');

fprintf('La distancia al Sol es: %d.', norm(x));

%Velocidad y su modulo

v=P(4)*(sqrt(mi))/(P(4)^(3/2)*(1-P(5)*cos(u0)))*[-sin(u0), sqrt(1-P(5)^2)*cos(u0)]
fprintf('\nEl modulo de la velocidad es: %d.', norm(v));

%Momento angular (en modulo)
%Dependiente del tiempo
c1=(P(4)^2)*sqrt(1-P(5)^2)*(sqrt(mi))/(P(4)^(3/2)*(1-P(5)*cos(u0)))*((cos(u0)-P(5))*cos(u0)+(sin(u0))^2);
%Usando las constantes
c2=sqrt(mi*P(4)*(1-P(5)^2));
fprintf('\nEl modulo del momento angular es %d calculado con la expresion dependiente del tiempo y\n %d calculado con la expresion que no depiende del tiempo.\n', c1, c2);

%Energia
%Dependiente del tiempo
h1=1/2*norm(v)^2-mi/norm(x);
%No dependiente dal tiempo
h2=mi^2*(P(5)^2-1)*(1/(2*(c2)^2));
fprintf('\nLa energia es %d calculada con la expresion dependiente del tiempo y\n %d calculada con la expresion que no depiende del tiempo.\n', h1, h2);

%Anomalia real dada la anomalia excentrica
exc=input('\nElegir un numero entre 0 y 2pi (la anomalia excentrica): '); 
costeta=(cos(exc)-P(5))/(1-P(5)*cos(exc));
sinteta=(sqrt(1-P(5)^2)*sin(exc))/(1-P(5)*cos(exc));
tgteta=costeta/sinteta;
teta=atan(tgteta);
fprintf('\nLa anomalia real obtenida en funcion de la anomalia excentrica es: %d.\n', teta);

%Bessel

sum = 0;

for n = 1:30
    bs = (2/n)*besselj(n,n*P(5))*sin(t*n*(sqrt(mi))/(P(4)^(3/2)));
    sum=sum+bs;
    n=n+1;
end

ubs=t*(sqrt(mi))/((P(4)^(3/2)))+sum
u0




