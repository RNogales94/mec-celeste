#include <iostream>
#include <cmath>
#include <math.h>
#include <string.h>

using namespace std;
//metodo de newtonraphson
// modificar vector y utilizar dos variables, una auxiliar, para que no te limite el hecho de 1000 posiciones del vector
double NR (double tolerancia,double epsilon, double periodo, double t, double tini) {
    double v[1000];
    double u0=0;
    double ji=(2*3.1415927*(t-tini))/periodo;
    int i=1;
    v[0]=u0;
    v[1]=u0-(u0-epsilon*sin(u0)-ji)/(1+epsilon*cos(u0));

    while(abs(v[i]-v[i-1])>tolerancia){
        v[i+1]=v[i]-((v[i]-epsilon*sin(v[i])-ji)/(1+epsilon*cos(v[i])));
        i++;
        }

     return (v[i]);

}
//runge kutta
double rungekutta (double a){
    double ar;




}
//pagina biblioteca funciones de bessel
/*https://msdn.microsoft.com/es-es/library/h7zkk1bz.aspx*/
double polinomiosdebessel(double tolerancia, double epsilon, double periodo, double t, double ti){
    double pi=3.141592;
    double sum=0;
    int i=1;
    while(i<1000){
        sum=sum+(2/i)*_jn(i,i*epsilon)*sin((2*pi*i*(t-ti))/periodo);
        i++;
    }
    return(((2*pi*(t-ti))/periodo)+sum);

}
//anomal�a real en funci�n de la exc�ntrica
double real(double u, double epsilon){
    double areal;
    areal=2*(atan(sqrt((1+epsilon)/(1-epsilon))*tan(u/2)));
    return(areal);
}

// calculo del vector a partir de la nomalia
double cv1(double u,double epsilon,double a){
    double i;
    i=a*(cos(u)-epsilon);
    return(i);

}
double cv2(double u,double epsilon, double a){
    double i;
    i=a*(sqrt(1-(epsilon*epsilon))*sin (u));
    return(i);

}
// calculo modulo vector
double cmv(double x, double y){
    double mod;
    mod= sqrt(x*x+y*y);
    return(mod);
}
// calculo derivada
double CD1(double tolerancia, double epsilon, double periodo, double t, double tini,double a){
    double i=t-1;
    double j;
    while ((t-i)>0.00002){
        j=(cv1(NR(tolerancia,epsilon,periodo,i,tini),epsilon,a)-cv1(NR(tolerancia,epsilon,periodo,t,tini),epsilon,a))/(i-t);
        i=i+0.001;
    }
    return(j);
}
double CD2(double tolerancia,double epsilon, double periodo, double t, double tini, double a){
    double i=t-1;
    double j=0;
    while ((t-i)>0.00002){
        j=(cv2(NR(tolerancia,epsilon,periodo,i,tini),epsilon,a)-cv2(NR(tolerancia,epsilon,periodo,t,tini),epsilon,a))/(i-t);
        i=i+0.001;
    }
    return(j);
}
double CDMA(double tolerancia, double epsilon, double periodo, double t, double tini,double a){
    double i=t-1;
    double j;
    while ((t-i)>0.00002){
        j=(NR(tolerancia,epsilon,periodo,i,tini)-NR(tolerancia,epsilon,periodo,t,tini))/(i-t);
        i=i+0.001;
    }
    return(j);
}
// calculo momento angular
double cma (double tol, double eps, double p, double t, double ti, double a){
    double angu;
    angu= cv1(NR (tol, eps, p, t, ti), eps, a)* CD2(tol, eps, p, t ,ti,a)-cv2(NR (tol, eps, p, t, ti), eps,a)* CD1(tol, eps, p, t ,ti,a);
    return(angu);
}
double calculomoduloangularu(double epsilon, double a, double u,double tolerancia, double periodo, double t, double tini){
    double ma;
    ma=a*a*sqrt(1-epsilon*epsilon)*abs(CDMA(tolerancia, epsilon, periodo, t,tini,a))*(1-epsilon*cos(u));
    return(ma);

}
double cma2 (double a, double epsilon, double mu){
    double ma;
    ma=sqrt(a*mu*((1-epsilon)*(1-epsilon)));
    return (ma);
}
// calculo energia total
double calculoenergiatotal (double tol,double eps, double p, double t, double ti, double M, double a,double mu){
	double et;
	et= (1/2)* (cmv(CD1(tol, eps, p, t ,ti,a),CD2(tol, eps, p, t,ti,a)))*(cmv(CD1(tol, eps, p, t ,ti,a),CD2(tol, eps, p, t,ti,a)))-(mu)/(cmv(cv1(NR(tol, eps, p, t,ti),eps,a),cv2(NR(tol,eps,p,t,ti),eps,a)));
    return(et);
}

//energia total en funcion de las ctes
double calculoenergiatotal2 (double a, double mu){
        double et;
        et=(-mu)/(a);
        return(et);
}



int main()
{
    double tolerancia=0.00000000001;
    double periodo;
    char planeta[20];
    double t;
    double epsilon;
    double tini=0;
    double M;
    double a;
    double mu;
    cout << "introduzca tiempo " << endl;
    cin >> t;
    cout << "introduzca planeta (todo en minuscula y sin tildes) " << endl;
    cin >> planeta;


    if (strcmp(planeta,"mercurio")==0){
    periodo=7600608;
    epsilon=0.206; // modulo excentricida/
    a=0.387; //
    mu=(4*3.141592*3.141592*a*a*a)/(periodo*periodo); // metros cubicos partido segundo cuadrados

    cout << "1.-la posicion del planeta es: " << cv1(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a) << " , " << cv2(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a) << endl;
    cout << "2.-la distancia al sol del planeta es: " << cmv(cv1(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a),cv2(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a)) << endl;
    cout << "3.-el vector velocidad es: " << CD1(tolerancia, epsilon, periodo, t ,tini,a) << " , "  << CD2(tolerancia, epsilon, periodo, t ,tini,a) << endl;
    cout << "3.-el modulo del vector velocidad es: " << cmv(CD1(tolerancia, epsilon, periodo, t ,tini,a),CD2(tolerancia, epsilon, periodo, t ,tini,a)) << endl;
    cout << "5.- la energia total usando la expresion del tiempo es: " << calculoenergiatotal(tolerancia, epsilon, periodo, t, tini,M,a, mu) << endl;
    cout << "5.- la energia total usando constantes es: " << calculoenergiatotal2(a,mu) << endl;
    cout << "6.- momento angular usando la expresion dependiente del tiempo: "<< cma(tolerancia, epsilon, periodo, t, tini,a) << endl;
    cout << "6.- momento angular usando la expresi�n dependiente del tiempo (en funcion de u): " << calculomoduloangularu(epsilon,a, NR (tolerancia, epsilon, periodo, t, tini), tolerancia,  periodo, t, tini) << endl;
    cout << "6.- momento angular mediante la constante calculada en funcion de los parametros del sistema:" << cma2(a,epsilon,mu) <<endl;
    cout << "7.- Anomalia real en funcion de la excentrica: " << real(NR (tolerancia, epsilon, periodo, t, tini),epsilon) << endl;
    cout << "8.- Anomalia excentrica bessel y comparacion con anomalia excentrica newton:" << endl;
    cout << "Newton raphson: " << NR (tolerancia, epsilon, periodo, t, tini)  << endl;
    cout << "Bessel: " << polinomiosdebessel(tolerancia, epsilon, periodo, t, tini) << endl;
    }
    if (strcmp(planeta,"venus")==0){
    periodo=14021280;
    epsilon=0.007; // modulo excentricida/
    a=0.723; //
    mu=(4*3.141592*3.141592*a*a*a)/(periodo*periodo); // metros cubicos partido segundo cuadrados

    cout << "1.-la posicion del planeta es: " << cv1(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a) << " , " << cv2(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a) << endl;
    cout << "2.-la distancia al sol del planeta es: " << cmv(cv1(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a),cv2(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a)) << endl;
    cout << "3.-el vector velocidad es: " << CD1(tolerancia, epsilon, periodo, t ,tini,a) << " , "  << CD2(tolerancia, epsilon, periodo, t ,tini,a) << endl;
    cout << "3.-el modulo del vector velocidad es: " << cmv(CD1(tolerancia, epsilon, periodo, t ,tini,a),CD2(tolerancia, epsilon, periodo, t ,tini,a)) << endl;
    cout << "5.- la energia total usando la expresion del tiempo es: " << calculoenergiatotal(tolerancia, epsilon, periodo, t, tini,M,a, mu) << endl;
    cout << "5.- la energia total usando constantes es: " << calculoenergiatotal2(a,mu) << endl;
    cout << "6.- momento angular usando la expresion dependiente del tiempo: "<< cma(tolerancia, epsilon, periodo, t, tini,a) << endl;
    cout << "6.- momento angular usando la expresi�n dependiente del tiempo (en funcion de u): " << calculomoduloangularu(epsilon,a, NR (tolerancia, epsilon, periodo, t, tini), tolerancia,  periodo, t, tini) << endl;
    cout << "6.- momento angular mediante la constante calculada en funcion de los parametros del sistema:" << cma2(a,epsilon,mu) <<endl;
    cout << "7.- Anomalia real en funcion de la excentrica: " << real(NR (tolerancia, epsilon, periodo, t, tini),epsilon) << endl;
    cout << "8.- Anomalia excentrica bessel y comparacion con anomalia excentrica newton:" << endl;
    cout << "Newton raphson: " << NR (tolerancia, epsilon, periodo, t, tini)  << endl;
    cout << "Bessel: " << polinomiosdebessel(tolerancia, epsilon, periodo, t, tini) << endl;
    }
    if (strcmp(planeta,"tierra")==0){
    periodo=31558464;
    epsilon=0.017; // modulo excentricida/
    a=1; //
    mu=(4*3.141592*3.141592*a*a*a)/(periodo*periodo); // metros cubicos partido segundo cuadrados

    cout << "1.-la posicion del planeta es: " << cv1(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a) << " , " << cv2(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a) << endl;
    cout << "2.-la distancia al sol del planeta es: " << cmv(cv1(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a),cv2(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a)) << endl;
    cout << "3.-el vector velocidad es: " << CD1(tolerancia, epsilon, periodo, t ,tini,a) << " , "  << CD2(tolerancia, epsilon, periodo, t ,tini,a) << endl;
    cout << "3.-el modulo del vector velocidad es: " << cmv(CD1(tolerancia, epsilon, periodo, t ,tini,a),CD2(tolerancia, epsilon, periodo, t ,tini,a)) << endl;
    cout << "5.- la energia total usando la expresion del tiempo es: " << calculoenergiatotal(tolerancia, epsilon, periodo, t, tini,M,a, mu) << endl;
    cout << "5.- la energia total usando constantes es: " << calculoenergiatotal2(a,mu) << endl;
    cout << "6.- momento angular usando la expresion dependiente del tiempo: "<< cma(tolerancia, epsilon, periodo, t, tini,a) << endl;
    cout << "6.- momento angular usando la expresi�n dependiente del tiempo (en funcion de u): " << calculomoduloangularu(epsilon,a, NR (tolerancia, epsilon, periodo, t, tini), tolerancia,  periodo, t, tini) << endl;
    cout << "6.- momento angular mediante la constante calculada en funcion de los parametros del sistema:" << cma2(a,epsilon,mu) <<endl;
    cout << "7.- Anomalia real en funcion de la excentrica: " << real(NR (tolerancia, epsilon, periodo, t, tini),epsilon) << endl;
    cout << "8.- Anomalia excentrica bessel y comparacion con anomalia excentrica newton:" << endl;
    cout << "Newton raphson: " << NR (tolerancia, epsilon, periodo, t, tini)  << endl;
    cout << "Bessel: " << polinomiosdebessel(tolerancia, epsilon, periodo, t, tini) << endl;
    }
    if (strcmp(planeta,"marte")==0){
    periodo=59355072;
    epsilon=0.093; // modulo excentricida/
    a=1.524; //
    mu=(4*3.141592*3.141592*a*a*a)/(periodo*periodo); // metros cubicos partido segundo cuadrados
    cout << "1.-la posicion del planeta es: " << cv1(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a) << " , " << cv2(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a) << endl;
    cout << "2.-la distancia al sol del planeta es: " << cmv(cv1(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a),cv2(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a)) << endl;
    cout << "3.-el vector velocidad es: " << CD1(tolerancia, epsilon, periodo, t ,tini,a) << " , "  << CD2(tolerancia, epsilon, periodo, t ,tini,a) << endl;
    cout << "3.-el modulo del vector velocidad es: " << cmv(CD1(tolerancia, epsilon, periodo, t ,tini,a),CD2(tolerancia, epsilon, periodo, t ,tini,a)) << endl;
    cout << "5.- la energia total usando la expresion del tiempo es: " << calculoenergiatotal(tolerancia, epsilon, periodo, t, tini,M,a, mu) << endl;
    cout << "5.- la energia total usando constantes es: " << calculoenergiatotal2(a,mu) << endl;
    cout << "6.- momento angular usando la expresion dependiente del tiempo: "<< cma(tolerancia, epsilon, periodo, t, tini,a) << endl;
    cout << "6.- momento angular usando la expresi�n dependiente del tiempo (en funcion de u): " << calculomoduloangularu(epsilon,a, NR (tolerancia, epsilon, periodo, t, tini), tolerancia,  periodo, t, tini) << endl;
    cout << "6.- momento angular mediante la constante calculada en funcion de los parametros del sistema:" << cma2(a,epsilon,mu) <<endl;
    cout << "7.- Anomalia real en funcion de la excentrica: " << real(NR (tolerancia, epsilon, periodo, t, tini),epsilon) << endl;
    cout << "8.- Anomalia excentrica bessel y comparacion con anomalia excentrica newton:" << endl;
    cout << "Newton raphson: " << NR (tolerancia, epsilon, periodo, t, tini)  << endl;
    cout << "Bessel: " << polinomiosdebessel(tolerancia, epsilon, periodo, t, tini) << endl;
    }
    if (strcmp(planeta,"jupiter")==0){
    periodo=374336640;
    epsilon=0.048; // modulo excentricida/
    a=5.203; //
    mu=(4*3.141592*3.141592*a*a*a)/(periodo*periodo); // metros cubicos partido segundo cuadrados

    cout << "1.-la posicion del planeta es: " << cv1(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a) << " , " << cv2(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a) << endl;
    cout << "2.-la distancia al sol del planeta es: " << cmv(cv1(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a),cv2(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a)) << endl;
    cout << "3.-el vector velocidad es: " << CD1(tolerancia, epsilon, periodo, t ,tini,a) << " , "  << CD2(tolerancia, epsilon, periodo, t ,tini,a) << endl;
    cout << "3.-el modulo del vector velocidad es: " << cmv(CD1(tolerancia, epsilon, periodo, t ,tini,a),CD2(tolerancia, epsilon, periodo, t ,tini,a)) << endl;
    cout << "5.- la energia total usando la expresion del tiempo es: " << calculoenergiatotal(tolerancia, epsilon, periodo, t, tini,M,a, mu) << endl;
    cout << "5.- la energia total usando constantes es: " << calculoenergiatotal2(a,mu) << endl;
    cout << "6.- momento angular usando la expresion dependiente del tiempo: "<< cma(tolerancia, epsilon, periodo, t, tini,a) << endl;
    cout << "6.- momento angular usando la expresi�n dependiente del tiempo (en funcion de u): " << calculomoduloangularu(epsilon,a, NR (tolerancia, epsilon, periodo, t, tini), tolerancia,  periodo, t, tini) << endl;
    cout << "6.- momento angular mediante la constante calculada en funcion de los parametros del sistema:" << cma2(a,epsilon,mu) <<endl;
    cout << "7.- Anomalia real en funcion de la excentrica: " << real(NR (tolerancia, epsilon, periodo, t, tini),epsilon) << endl;
    cout << "8.- Anomalia excentrica bessel y comparacion con anomalia excentrica newton:" << endl;
    cout << "Newton raphson: " << NR (tolerancia, epsilon, periodo, t, tini)  << endl;
    cout << "Bessel: " << polinomiosdebessel(tolerancia, epsilon, periodo, t, tini) << endl;
    }
    if (strcmp(planeta,"saturno")==0){
    periodo=929577600;
    epsilon=0.056; // modulo excentricida/
    a=9.546; //
    mu=(4*3.141592*3.141592*a*a*a)/(periodo*periodo); // metros cubicos partido segundo cuadrados

    cout << "1.-la posicion del planeta es: " << cv1(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a) << " , " << cv2(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a) << endl;
    cout << "2.-la distancia al sol del planeta es: " << cmv(cv1(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a),cv2(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a)) << endl;
    cout << "3.-el vector velocidad es: " << CD1(tolerancia, epsilon, periodo, t ,tini,a) << " , "  << CD2(tolerancia, epsilon, periodo, t ,tini,a) << endl;
    cout << "3.-el modulo del vector velocidad es: " << cmv(CD1(tolerancia, epsilon, periodo, t ,tini,a),CD2(tolerancia, epsilon, periodo, t ,tini,a)) << endl;
    cout << "5.- la energia total usando la expresion del tiempo es: " << calculoenergiatotal(tolerancia, epsilon, periodo, t, tini,M,a, mu) << endl;
    cout << "5.- la energia total usando constantes es: " << calculoenergiatotal2(a,mu) << endl;
    cout << "6.- momento angular usando la expresion dependiente del tiempo: "<< cma(tolerancia, epsilon, periodo, t, tini,a) << endl;
    cout << "6.- momento angular usando la expresi�n dependiente del tiempo (en funcion de u): " << calculomoduloangularu(epsilon,a, NR (tolerancia, epsilon, periodo, t, tini), tolerancia,  periodo, t, tini) << endl;
    cout << "6.- momento angular mediante la constante calculada en funcion de los parametros del sistema:" << cma2(a,epsilon,mu) <<endl;
    cout << "7.- Anomalia real en funcion de la excentrica: " << real(NR (tolerancia, epsilon, periodo, t, tini),epsilon) << endl;
    cout << "8.- Anomalia excentrica bessel y comparacion con anomalia excentrica newton:" << endl;
    cout << "Newton raphson: " << NR (tolerancia, epsilon, periodo, t, tini)  << endl;
    cout << "Bessel: " << polinomiosdebessel(tolerancia, epsilon, periodo, t, tini) << endl;
    }


    if (strcmp(planeta,"urano")==0){
    periodo=2651356800;
    epsilon=0.047; // modulo excentricida/
    a=19.20; //
    mu=(4*3.141592*3.141592*a*a*a)/(periodo*periodo); // metros cubicos partido segundo cuadrados

    cout << "1.-la posicion del planeta es: " << cv1(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a) << " , " << cv2(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a) << endl;
    cout << "2.-la distancia al sol del planeta es: " << cmv(cv1(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a),cv2(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a)) << endl;
    cout << "3.-el vector velocidad es: " << CD1(tolerancia, epsilon, periodo, t ,tini,a) << " , "  << CD2(tolerancia, epsilon, periodo, t ,tini,a) << endl;
    cout << "3.-el modulo del vector velocidad es: " << cmv(CD1(tolerancia, epsilon, periodo, t ,tini,a),CD2(tolerancia, epsilon, periodo, t ,tini,a)) << endl;
    cout << "5.- la energia total usando la expresion del tiempo es: " << calculoenergiatotal(tolerancia, epsilon, periodo, t, tini,M,a, mu) << endl;
    cout << "5.- la energia total usando constantes es: " << calculoenergiatotal2(a,mu) << endl;
    cout << "6.- momento angular usando la expresion dependiente del tiempo: "<< cma(tolerancia, epsilon, periodo, t, tini,a) << endl;
    cout << "6.- momento angular usando la expresi�n dependiente del tiempo (en funcion de u): " << calculomoduloangularu(epsilon,a, NR (tolerancia, epsilon, periodo, t, tini), tolerancia,  periodo, t, tini) << endl;
    cout << "6.- momento angular mediante la constante calculada en funcion de los parametros del sistema:" << cma2(a,epsilon,mu) <<endl;
    cout << "7.- Anomalia real en funcion de la excentrica: " << real(NR (tolerancia, epsilon, periodo, t, tini),epsilon) << endl;
    cout << "8.- Anomalia excentrica bessel y comparacion con anomalia excentrica newton:" << endl;
    cout << "Newton raphson: " << NR (tolerancia, epsilon, periodo, t, tini)  << endl;
    cout << "Bessel: " << polinomiosdebessel(tolerancia, epsilon, periodo, t, tini) << endl;
    }
    if (strcmp(planeta,"neptuno")==0){


    periodo=5251737600;
    epsilon=0.009; // modulo excentricida/
    a=30.09; //
    mu=(4*3.141592*3.141592*a*a*a)/(periodo*periodo); // metros cubicos partido segundo cuadrados

    cout << "1.-la posicion del planeta es: " << cv1(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a) << " , " << cv2(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a) << endl;
    cout << "2.-la distancia al sol del planeta es: " << cmv(cv1(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a),cv2(NR (tolerancia, epsilon, periodo, t, tini), epsilon,a)) << endl;
    cout << "3.-el vector velocidad es: " << CD1(tolerancia, epsilon, periodo, t ,tini,a) << " , "  << CD2(tolerancia, epsilon, periodo, t ,tini,a) << endl;
    cout << "3.-el modulo del vector velocidad es: " << cmv(CD1(tolerancia, epsilon, periodo, t ,tini,a),CD2(tolerancia, epsilon, periodo, t ,tini,a)) << endl;
    cout << "5.- la energia total usando la expresion del tiempo es: " << calculoenergiatotal(tolerancia, epsilon, periodo, t, tini,M,a, mu) << endl;
    cout << "5.- la energia total usando constantes es: " << calculoenergiatotal2(a,mu) << endl;
    cout << "6.- momento angular usando la expresion dependiente del tiempo: "<< cma(tolerancia, epsilon, periodo, t, tini,a) << endl;
    cout << "6.- momento angular usando la expresi�n dependiente del tiempo (en funcion de u): " << calculomoduloangularu(epsilon,a, NR (tolerancia, epsilon, periodo, t, tini), tolerancia,  periodo, t, tini) << endl;
    cout << "6.- momento angular mediante la constante calculada en funcion de los parametros del sistema:" << cma2(a,epsilon,mu) <<endl;
    cout << "7.- Anomalia real en funcion de la excentrica: " << real(NR (tolerancia, epsilon, periodo, t, tini),epsilon) << endl;
    cout << "8.- Anomalia excentrica bessel y comparacion con anomalia excentrica newton:" << endl;
    cout << "Newton raphson: " << NR (tolerancia, epsilon, periodo, t, tini)  << endl;
    cout << "Bessel: " << polinomiosdebessel(tolerancia, epsilon, periodo, t, tini) << endl;
    }
}
