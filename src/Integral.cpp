#include "Utilities.hpp"
#include <math.h>
#include "Integral.hpp"

const double EPSq = 1.0e-16;
const int JMAX    = 40;
const int JMAXP   = (JMAX+1);
const int K       = 12;

/****** Função de interpolação necessária para qromb ******/
void polint(double xa[], double ya[], int n, double x, double *y, double *dy) {
   
   /* Declaração das variáveis */
   int i,m,ns=1;
   double den,dif,dift,ho,hp,w;
   double *c,*d;
   
   dif=fabs(x-xa[1]);
   c=vector<double>(1,n);
   d=vector<double>(1,n);
   for (i=1;i<=n;i++) {
      if ( (dift=fabs(x-xa[i])) < dif) {	
	 ns=i;
	 dif=dift;
      }
      c[i]=ya[i];
      d[i]=ya[i];
   }
   *y=ya[ns--];
   for (m=1;m<n;m++) {
      for (i=1;i<=n-m;i++) {
	 ho=xa[i]-x;
	 hp=xa[i+m]-x;
	 w=c[i+1]-d[i];
	 if ( (den=ho-hp) == 0.0) error("Error in routine polint");
	 den=w/den;
	 d[i]=hp*den;
	 c[i]=ho*den;
      }
      *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
   }
   
   free_vector(d,1,n);
   free_vector(c,1,n);
}


  
/****** Função necessária para a integração qromb2 ******/
// Used for functions that receive Cosmological parameters through 'Cosmology' object. trapzd2.
double trapzd(double (*func)(param,double), double a, double b, int n, param p) {
   
   /* Declaração das variáveis */
   double x,tnm,sum,del;
   static double s;
   int it,j;
   
   if (n == 1) {
     return (s=0.5*(b-a)*((*func)(p,a)+(*func)(p,b)));
   }
   else {
      for (it=1,j=1;j<n-1;j++) it <<= 1;
      tnm=it;
      del=(b-a)/tnm;
      x=a+0.5*del;
      for (sum=0.0,j=1;j<=it;j++,x+=del) sum += (*func)(p,x);
      s=0.5*(s+(b-a)*sum/tnm);
      return s;
   }
}

/****** Método de integracao de Romberg para funções boazinhas ******/
// This version works with functions that need cosmological parameters. qromb2.
double qromb(double (*func)(param,double), double a, double b, param p) {
   
   /* Declaração de funçõe utilizadas */
   void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
   double trapzd(double (*func)(param , double), double a, double b, int n, param p);
   /* Declaração das variáveis utilizadas */
   double ss,dss, EPSqi;
   double s[JMAXP],h[JMAXP+1];
   int j, Ki, JMAXi;
   
   Ki=K;
   EPSqi=EPSq;
   JMAXi=JMAX;
      
   h[1]=1.0;
   for (j=1; j<=JMAXi; j++) {
     s[j]=trapzd(func,a,b,j,p);
      if (j >= Ki) {
	 polint(&h[j-Ki],&s[j-Ki],Ki,0.0,&ss,&dss);
	 if (fabs(dss) <= EPSqi*fabs(ss)) return ss;
      }
      h[j+1]=0.25*h[j];
   }
   error("Too many steps in routine qromb2");
   
   return 0.0;
}

/****** Função necessária para a integração qromb ******/
// For qromb3, used for integrating functions with extra 'double' parameter 'p'. trapzd3.
double trapzd(double (*func)(double, double), double a, double b, int n, double p) {
   
   /* Declaração das variáveis */
   double x,tnm,sum,del;
   static double s;
   int it,j;
   
   if (n == 1) {
     return (s=0.5*(b-a)*((*func)(p,a)+(*func)(p,b)));
   }
   else {
      for (it=1,j=1;j<n-1;j++) it <<= 1;
      tnm=it;
      del=(b-a)/tnm;
      x=a+0.5*del;
      for (sum=0.0,j=1;j<=it;j++,x+=del) sum += (*func)(p,x);
      s=0.5*(s+(b-a)*sum/tnm);
      return s;
   }
}

/****** Método de integracao de Romberg para funções boazinhas ******/
// For functions with an extra 'double' parameter 'p'. qromb3.
double qromb(double (*func)(double, double), double a, double b, double p) {
   /* Declaração de funçõe utilizadas */
   void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
   double trapzd(double (*func)(double, double), double a, double b, int n, double p);
   /* Declaração das variáveis utilizadas */
   double ss,dss, EPSqi;
   double s[JMAXP],h[JMAXP+1];
   int j, Ki, JMAXi;
   
   Ki=K;
   EPSqi=EPSq;
   JMAXi=JMAX;
      
   h[1]=1.0;
   for (j=1; j<=JMAXi; j++) {
     s[j]=trapzd(func,a,b,j,p);
      if (j >= Ki) {
	 polint(&h[j-Ki],&s[j-Ki],Ki,0.0,&ss,&dss);
	 if (fabs(dss) <= EPSqi*fabs(ss)) return ss;
      }
      h[j+1]=0.25*h[j];
   }
   error("Too many steps in routine qromb3");
   
   return 0.0;
}


/****** Função necessária para a integração qromb ******/
// For functions with extra double AND parameters.
double trapzd(double (*func)(double, double, param), double a, double b, int n, double z0, param p) {
   
   /* Declaração das variáveis */
   double x,tnm,sum,del;
   static double s;
   int it,j;
   
   if (n == 1) {
     return (s=0.5*(b-a)*((*func)(a,z0,p)+(*func)(b,z0,p)));
   }
   else {
      for (it=1,j=1;j<n-1;j++) it <<= 1;
      tnm=it;
      del=(b-a)/tnm;
      x=a+0.5*del;
      for (sum=0.0,j=1;j<=it;j++,x+=del) sum += (*func)(x,z0,p);
      s=0.5*(s+(b-a)*sum/tnm);
      return s;
   }
}

// Almost just a copy of above to allow functions with parameters as first entries to be integrated:
double trapzd5(double (*func)(param, double, double), double a, double b, int n, double z0, param p) {
   
   /* Declaração das variáveis */
   double x,tnm,sum,del;
   static double s;
   int it,j;
   
   if (n == 1) {
     return (s=0.5*(b-a)*((*func)(p,a,z0)+(*func)(p,b,z0)));
   }
   else {
      for (it=1,j=1;j<n-1;j++) it <<= 1;
      tnm=it;
      del=(b-a)/tnm;
      x=a+0.5*del;
      for (sum=0.0,j=1;j<=it;j++,x+=del) sum += (*func)(p,x,z0);
      s=0.5*(s+(b-a)*sum/tnm);
      return s;
   }
}

/****** Método de integracao de Romberg para funções boazinhas ******/
// For functions with extra double AND parameters.
double qromb(double (*func)(double, double, param), double a, double b, double z0, param p) {
   /* Declaração de funçõe utilizadas */
   void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
   double trapzd(double (*func)(double, double, param), double a, double b, int n, double z0, param p);
   /* Declaração das variáveis utilizadas */
   double ss,dss, EPSqi;
   double s[JMAXP],h[JMAXP+1];
   int j, Ki, JMAXi;
   
   Ki=K;
   EPSqi=EPSq;
   JMAXi=JMAX;
      
   h[1]=1.0;
   for (j=1; j<=JMAXi; j++) {
     s[j]=trapzd(func,a,b,j,z0,p);
      if (j >= Ki) {
	 polint(&h[j-Ki],&s[j-Ki],Ki,0.0,&ss,&dss);
	 if (fabs(dss) <= EPSqi*fabs(ss)) return ss;
      }
      h[j+1]=0.25*h[j];
   }
   error("Too many steps in routine qromb4");
   
   return 0.0;
}

// Just a copy of above to allow functions with parameters as first entries to be integrated:
double qromb5(double (*func)(param, double, double), double a, double b, double z0, param p) {
  /* Declaração de funçõe utilizadas */
   void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
   double trapzd(double (*func)(double, double, param), double a, double b, int n, double z0, param p);
   /* Declaração das variáveis utilizadas */
   double ss,dss, EPSqi;
   double s[JMAXP],h[JMAXP+1];
   int j, Ki, JMAXi;
   
   Ki=K;
   EPSqi=EPSq;
   JMAXi=JMAX;
      
   h[1]=1.0;
   for (j=1; j<=JMAXi; j++) {
     s[j]=trapzd5(func,a,b,j,z0,p);
      if (j >= Ki) {
	 polint(&h[j-Ki],&s[j-Ki],Ki,0.0,&ss,&dss);
	 if (fabs(dss) <= EPSqi*fabs(ss)) return ss;
      }
      h[j+1]=0.25*h[j];
   }
   error("Too many steps in routine qromb5");
   
   return 0.0;
}
