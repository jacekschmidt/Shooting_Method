#include <iostream>
#include "mvector.h"
#include <cmath>
#include <fstream>

//base class
class MFunction
{
    public:
        virtual MVector operator()(const double& x, const MVector& y) = 0;
};

//derived f1 class
class FunctionF1 : public MFunction
{
    public:
        virtual MVector operator()(const double& x, const MVector& y)
        {
        MVector temp(2);
        temp[0] = y[0] + x*y[1];
        temp[1] = x*y[0] - y[1];
        return temp;
        }
};

//derived class F2
class FunctionF2 : public MFunction
{
    public:
        virtual MVector operator()(const double& x, const MVector& y)
        {
        MVector temp(2);
        temp[0] = x;
        temp[1] =y[1];
        return temp;
        }
};

//derived class F3
class FunctionF3 : public MFunction
{
    public:
        virtual MVector operator()(const double& x, const MVector& y)
        {
        MVector temp(1);
        temp[0] = 2*x;
        return temp;
        }
};

// derived eqn1p5 class
class Eqn1p5Derivs : public MFunction
{
    public:
        // constructor to initialise kappa
        Eqn1p5Derivs() {kappa=1.0;}
        MVector operator()(const double& x,const MVector& y)
        {
        MVector temp(4);
        temp[0] = y[1];
        temp[1] = -kappa*y[1] - x*y[0];
        temp[2] = y[3];
        temp[3] = -kappa*y[3] - x*y[2];
        return temp;
        }
        void SetKappa(double k) {kappa=k;} // change kappa
    private:
        double kappa; // class member variable, accessible within
        // all Eqn1p5Derivs member functions
};

// derived Falkner Skan class
class FalknerSkan : public MFunction
{
    public:
        // constructor to initialise kappa
        FalknerSkan() {beta=0.0;}
        MVector operator()(const double& x,const MVector& y)
        {
        MVector temp(6);
        temp[0] = y[1];
        temp[1] = y[2];
        temp[2] = beta*(y[1]*y[1]-1)-y[0]*y[2];
        temp[3] = y[4];
        temp[4] = y[5];
        temp[5] = 2*beta*y[1]*y[4]-y[2]*y[3]-y[0]*y[5];
        return temp;
        }
        void SetBeta(double b) {beta=b;} // change kappa
    private:
        double beta; // class member variable, accessible within
        // all Eqn1p5Derivs member functions
};

// derived quadratic tester class
class quadratictester : public MFunction
{
    public:
        MVector operator()(const double& x,const MVector& y)
        {
        MVector temp(9);
        temp[0] = y[1];
        temp[1] = y[2];
        temp[2] = 0;
        temp[3] = y[4];
        temp[4] = y[5];
        temp[5] = 0;
        temp[6] = y[7];
        temp[7] = y[8];
        temp[8] = 0;
        return temp;
        }
};




// Declaration for an Euler scheme ODE solver
int EulerSolve(int steps, double a, double b, MVector &Y, MFunction &F)
{
    std::ofstream demofile;
    demofile.open("Euler_ODE_solve");
    if (!demofile) return 1;
    double h=(b-a)/steps;
    if (b<a) return 1;
    demofile<<a<<" , "<<Y<<"\n";
    for(int i=0; i<steps;i++)
    {
        Y=Y+h*F(a,Y);
        a+=h;
        demofile<<a<<" , "<<Y<<"\n";
    }
    demofile.close();
    return 0;
}


// Runge Kutta Solver
int RungeKuttaSolve(int steps, double a, double b, MVector &Y, MFunction &F)
{   
    std::ofstream demofile;
    demofile.open("Runge_Kutta_ODE_solve");
    if (!demofile) return 1;
    MVector k1 ,k2, k3, k4;
    double h=(b-a)/steps;
    if (b<a) return 1;
    demofile<<a<<" , "<<Y<<"\n";
    for(int i=0; i<steps;i++)
    {
        k1=F(a,Y);
        k2=F(a+(h/2),Y+(h/2)*k1);
        k3=F(a+(h/2),Y+(h/2)*k2);
        k4=F(a+h,Y+h*k3);
        Y=Y+(h/6)*(k1+2*k2+2*k3+k4);
        a+=h;
        demofile<<a<<" , "<<Y<<"\n";
        std::cout<<Y[0]<<std::endl;
    }
    demofile.close();
    return 0;
}

MVector FalknerSkanGuessr(double beta, double guess, int end, int maxNewtonSteps)
{
    FalknerSkan f;
    f.SetBeta(beta);
    MVector y;
    int maxnumintstep=end*100;
    double tol=1e-8, phi;     
    for (int i=0; i<maxNewtonSteps;i++)
    {
        y={0,0,guess,0,0,1};
        RungeKuttaSolve(maxnumintstep,0,end,y,f);
        phi=y[1]-1;
        double phidash=y[4];
        if (std::abs(phi)<tol) break;
        guess -= phi/phidash;
    }
    MVector v={beta,guess};
    return v;
}

MVector quadraticguessr(double guessp, double guessq, int descentSteps, double nu, double tol)
{
    quadratictester f;
    MVector y;
    double phi=1.0, phi1, phi2, phi1p, phi1q, phi2p, phi2q;
    double j=0;
    for(int i=0;i<descentSteps;i++)
    {
        phi=0.0;
        y={1,guessp,guessq,0,1,0,0,0,1};
        RungeKuttaSolve(0.5*100,0,0.5,y,f);
        phi1=y[0];
        phi+=phi1*phi1;
        phi1p=y[3];
        phi1q=y[6];
        y={1,guessp,guessq,0,1,0,0,0,1};
        RungeKuttaSolve(100,0,1,y,f);
        phi2=(y[0]-2);
        phi+=phi2*phi2;
        phi2p=y[3];
        phi2q=y[6];
        if (std::abs(phi)<tol) break;
        guessp-=nu*2*(phi1*phi1p+phi2*phi2p);
        guessq-=nu*2*(phi1*phi1q+phi2*phi2q);
        j++;
    }
    MVector V={guessp,guessq,phi,j};
    return V;
}

int main()
{   /*report q1  
    {
        FunctionF2 f;
        MVector y={0,1};
        EulerSolve(100,0,1,y,f);
        y={0,1};
        RungeKuttaSolve(100,0,1,y,f);
        std::cout<<FalknerSkanGuessr(0.1, 0.4, 5, 5)<<std::endl;
    }
    */
    /*report  q2
    {
        FunctionF2 f;
        MVector y={0,1};
        for(int i=1;i<=100;i++)
        {
            RungeKuttaSolve(i,0,1,y,f);
            std::cout.precision(12);
            std::cout<<y[1]<<std::endl;
            y={0,1};
        }
    }
    */
    /* later
   {
    double g=0.47, beta=0;
    double h=1/100.0;
    
    MVector VIn={5,10,15,20,25,30,50,75,100,125,150,160,170,180,190,200};
    for(int i=0;i<=100;i++)
    {
        for(int j=0;j<VIn.size();j++)
        {
            g=FalknerSkanGuessr(beta,g,VIn[j],10)[1];
        }
        if (g>10) break;
        std::cout.precision(18);
        std::cout<< FalknerSkanGuessr(beta,g,200,5)[1]<<std::endl;
        beta+=h;
    }
   }
    */
   /*residuals
   {
    FalknerSkan f;
    f.SetBeta(10.0/11.0);
    MVector y;
    int end=20;
    int maxnumintstep=end*100;
    double guess=1;
    double h=0.001;
    for ( int i=0; i<=500;i++)
    {
        y={0,0,guess,0,0,1};
    RungeKuttaSolve(maxnumintstep,0,end,y,f);
    if(std::abs(y[1]-1)<100)
    {
        std::cout<<y[1]-1<<std::endl;
    }
    guess+=h;
    }
   }
    */
   /*f''
   {
    double guess0=0.4695999884;
    double guessthird=0.8021255928;
    double guesselevth=1.182817236074818;
    FalknerSkan f;
    f.SetBeta(1.0/3.0);
    MVector y={0,0,guessthird+0.25*guessthird,0,0,1};
    int end=5;
    int maxnumintstep=end*100;
    RungeKuttaSolve(maxnumintstep,0,end,y,f);
   }
    */
   /*extra
    {
        quadratictester f;
        MVector y;
        double phi=0.0, phi1, phi2, phi1p, phi1q, phi2p, phi2q, end1=0.5, end2=1.0, guessp=1, guessq=2, nu=0.01;
        y={1,guessp,guessq,0,1,0,0,0,1};
        RungeKuttaSolve(end1*100,0,end1,y,f);
        phi1=y[0];
        phi+=phi1*phi1;
        phi1p=y[3];
        phi1q=y[6];
        y={1,guessp,guessq,0,1,0,0,0,1};
        RungeKuttaSolve(end2*100,0,end2,y,f);
        phi2=(y[0]-2);
        phi+=phi2*phi2;
        phi2p=y[3];
        phi2q=y[6];
        std::cout<<guessp<<" , "<<guessq<<" , "<<phi<<std::endl;
        guessp-=nu*2*(phi1*phi1p+phi2*phi2p);
        guessq-=nu*2*(phi1*phi1q+phi2*phi2q);
        phi=0;
        y={1,guessp,guessq,0,1,0,0,0,1};
        RungeKuttaSolve(end1*100,0,end1,y,f);
        phi1=y[0];
        phi+=phi1*phi1;
        phi1p=y[3];
        phi1q=y[6];
        y={1,guessp,guessq,0,1,0,0,0,1};
        RungeKuttaSolve(end2*100,0,end2,y,f);
        phi2=(y[0]-2);
        phi+=phi2*phi2;
        phi2p=y[3];
        phi2q=y[6];
        std::cout<<guessp<<" , "<<guessq<<" , "<<phi<<std::endl;

    }
1,3,7*/
{
    //std::cout<<quadraticguessr(1,1,700,0.1,0.001)<<std::endl;
    quadratictester f;
    MVector y={1,-4.20448,10.2469,0,1,0,0,0,1};
    RungeKuttaSolve(100,0,1,y,f);
}

    return 0;
}