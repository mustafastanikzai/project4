#include <iostream>
#include <math.h>
#include <armadillo>
#include <fstream>

using namespace arma;
using namespace std;

double g(double x);

int main()
{
    ofstream myfile;
    ofstream myfile2;
    ofstream myfile3;

    ofstream tilplot;
    tilplot.open("Data.txt");

    myfile.open("ForwardEuler.txt");
    myfile2.open("BackwardEuler.txt");
    myfile3.open("Crank_Nicolson.txt");

    // Deklarasjoner
    int N;
    int n; double alpha, d, D, dx, dt, T;
    double t;
    int timesteps;

    // Input
    timesteps = 30000; // = t2
    //timesteps = 30; // = t1

    n = 15; // antall grid punkter
    d = 1; // avstanden
    D = 1;
    N = 10000;

    // Kalkulasjoner
    dx = (double) d/(n-1);
    dt = dx*dx/50;

    T = dt * timesteps; // antall sekunder

    alpha = dt/(dx*dx);

    mat A = zeros(n,n); //matrise A med størrelse nxn
    vec Vnew = zeros(n);
    vec Vold = zeros(n);

    A(0,0) = 1-2*alpha;
    A(0,1) = alpha;

    for (int i=1; i < n-1; i++)
    {
        A(i,i) = 1-2*alpha;
        A(i,i+1) = alpha;
        A(i,i-1) = alpha;
    }
    A(n-1,n-1) = 1-2*alpha;
    A(n-1,n-2) = alpha;

    // Initialiserer u(x=0,t)
    Vold(0) = 1;
    Vold(n-1) = 0;
    Vnew(0) = 1;
    Vnew(n-1) = 0;

    ////////////////////////////////////////////////////////////
    // Forward Euler
    t=0;
    while (t < T)
    {
        for (int i=0; i<n; i++)
        {
            myfile << Vold(i) << " ";
        }

        //Vnew(0) = alpha * Vold(1) + (1-2*alpha)*Vold(0);

        for (int i = 1; i<n-1; i++)
        {
            Vnew(i) = (alpha*Vold(i+1) + (1-2*alpha) * Vold(i) + alpha*Vold(i-1));
        }
        t += dt;
        Vold = Vnew;
        myfile << endl;
    }

    vec Forward_Euler_Vec = Vnew;

    //////////////////////////////////////////////////////////
    // Backward Euler
    for (int i=0; i<n; i++)
    {
        Vold(i) = 0;
        Vnew(i) = 0;
    }

    // Vold er nå høyre side, Vnew venstre side av v_new = A^-1 v_old

    // Initialiser systemet
    Vold(0) = 1.0;
    Vnew(0) = 1.0;
    vec gamma = zeros(n);
    vec Valt = zeros(n);

    double bet;

    t = 0;
    while (t < T)
    {
        for (int i=0; i<n; i++)
        {
            myfile2 << Vold(i) << " ";
        }

        bet = 1 + 2 * alpha;
        for (int i = 1; i < n; i++)
        {
            gamma(i) = -alpha/bet; //(1+ 2*alpha);
            bet = (1 + 2*alpha) + alpha*gamma(i);
            if (bet == 0.0){cout << ":(" << endl;}
            Vnew(i) = (Vold(i) + alpha*Vnew(i-1))/bet;
        }

        for (int i=n-2; i > 0; i--)
        {
            Vnew(i) -= gamma(i+1)*Vnew(i+1);
        }
        Vnew(n-1) = 0;
        Vold = Vnew;
        myfile2 << endl;
        t += dt;
    }

    vec Backward_Euler_Vec = Vnew;

    ////////////////////////////////////////////////////////////
    // Crank Nicolson
    for (int i=0; i<n; i++)
    {
        Vnew(i) = 0;
        Vold(i) = 0;
        gamma(i) = 0;
    }

    Vold(0) = 1;
    Vnew(0) = 1;

    t = 0;
    vec r = zeros(n);

    while (t < T)
    {
        r(0) = Vold(0)*(2-2*alpha) + alpha*Vold(1);
        for (int i = 1; i< n-1; i++)
        {
            r(i) = Vold(i-1)*alpha + (2 - 2 * alpha)*Vold(i) + alpha * Vold(i+1);
        }
        r(n-1) = Vold(n-2)*alpha + (2-2*alpha)*Vold(n-1);

        // Så tridiagonal biten
        bet = 2+2*alpha;

        for (int i = 1; i < n; i++)
        {
            gamma(i) = -alpha/bet;
            bet = (2+2*alpha) + alpha*gamma(i);
            if (bet == 0.0){cout << ":(" << endl;}
            Vnew(i) = (r(i) + alpha * Vnew(i-1))/bet;
        }

        for (int i = n-2; i > 0; i--)
        {
            Vnew(i) -= gamma(i+1)*Vnew(i+1);
        }
        Vnew(n-1) = 0;

        Vold = Vnew;
        t += dt;
    }

    vec Crank_Nicolson_Vec = Vnew;

    //////////////////////////////////////////////////////
    // Eksakt
    for (int i = 0; i < n; i++)
    {
        Vnew(i) = 0;
        Vold(i) = 0;
    }

    t = T - dt; // Tar eksakt funksjon i samme tidssteg\

    for (int i = 0; i< n; i++)
    {
       //


        double s = 0.0;
        for (int j = 1; j <= N; j++)
        {
            s += sin(j*acos(-1.0)*i*dx) / j * exp(-j*j*acos(-1.0)*acos(-1.0)*t);
        }
        Vnew(i) = 1 - i*dx - 2.0/acos(-1.0) * s;

    }
    Vnew(n-1) = 0;

    vec Eksakt_Vec = Vnew;

    // Sammenligne resultater

    for (int i=0; i<n; i++)
    {
        tilplot << Backward_Euler_Vec(i) << " "
             << Forward_Euler_Vec(i) << " "
             << Crank_Nicolson_Vec(i) << " "
             << Eksakt_Vec(i) << endl;
        cout << Backward_Euler_Vec(i) << " " << Forward_Euler_Vec(i) << " " << Crank_Nicolson_Vec(i) << " " << Eksakt_Vec(i) << endl;
    }

    cout << "Data printes til Data.txt" << endl;
    cout << "Alpha brukt: " << alpha << ", alpha bør være under 0.5 helst en del under." << endl;

    myfile.close();
    myfile2.close();
    myfile3.close();
    tilplot.close();

    // vec Truncation_Error = sum(abs(Eksakt_Vec - Forward_Euler_Vec));

    return 0;
}

double g(double x)
{
    return 1; // :D
}
