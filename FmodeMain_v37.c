#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#define pi 3.14159265358979
#define G 6.673e-11
#define c 2.99792458e8
#define M 1.989e30
#define Cri 1e-7
#define Gc2 7.42471382405e-28
#define Gc4 8.26110825251e-45

static double a[3000],b[3000],l=2;
static int j,scr;
char dysm[6]="/-\\/-\\";
FILE *inf,*outf;
double temp_rhoc,temp_M,temp_R,temp_freq,temp_dpt;
double E_core,Emax,Emin;
double *mfile,*pfile,*rhofile,*Bfile,*Bcor,*Wfile,*Wfile1,*Wfile2,*Vfile,*Vfile1,*Vfile2;
struct EoS
	{
		char Name[36],tName[30];
	};
struct EoS list[50001];

double che(e)
double e;
{
double le;
double p,lp;
int i=0;
le=log10(e);
i=j-1;
for(;i>0;i--)
{
if(le>a[i-1])
  {
  lp=b[i-1]*(le-a[i])/(a[i-1]-a[i])+b[i]*(le-a[i-1])/(a[i]-a[i-1]);
  p=pow(10,lp);
  return(p);
  }
}
return(0);
}

double ch(p)
double p;
{
double lp;
double e,le;
int i=0;
lp=log10(p);
i=j-1;
for(;i>0;i--)
{
if(lp>b[i-1])
  {
  le=a[i-1]*(lp-b[i])/(b[i-1]-b[i])+a[i]*(lp-b[i-1])/(b[i]-b[i-1]);
  e=pow(10,le);
  return(e);
  }
}
return(0);
}

double fp(r,p,e,m)
double r,p,e,m;
{
return(-G*(e+p/c/c)*(m+4*pi*r*r*r*p/c/c)/(r*r-2*r*m*Gc2));
}

double fm(r,e)
double r,e;
{
return(4*pi*r*r*e);
}

double Bf(r,p,m,B)
double r,p,m,B;
{
return(2*Gc2/r/r*(m+4*pi*r*r*r*p/c/c)*B/(1-2*Gc2*m/r));
}

double DH1(r,m,A,p,e,H1,H0,K,V)
double r,m,A,p,e,H1,H0,K,V;
{
return(-1/r*(l+1+2*m*A/r+4*pi*r*r*A*(p-e))*H1+1/r*A*(H0+K-16*pi*(e+p)*V));
}

double DK(r,H0,H1,Dv,K,e,p,A,W)
double r,H0,H1,Dv,K,e,p,A,W;
{
return(1/r*H0+0.5*l*(l+1)/r*H1-((l+1)/r-0.5*Dv)*K-8*pi*(e+p)*sqrt(A)/r*W);
}

double DW(r,W,A,gamma,p,B,X,V,H0,K)
double r,W,A,gamma,p,B,X,V,H0,K;
{
return(-(l+1)/r*W+r*sqrt(A)*(1/gamma/p/sqrt(B)*X-l*(l+1)/r/r*V+0.5*H0+K));
}

double DX(r,X,e,p,B,Dv,H0,w,H1,K,V,A,F,W)
double r,X,e,p,B,Dv,H0,w,H1,K,V,A,F,W;
{
return(-l/r*X+(e+p)*sqrt(B)*(0.5*(1/r-0.5*Dv)*H0+0.5*(r*w*w/B+0.5*l*(l+1)/r)*H1+0.5*(1.5*Dv-1/r)*K-0.5*l*(l+1)*Dv/r/r*V-1/r*(4*pi*(e+p)*sqrt(A)+w*w*sqrt(A)/B-0.5*r*r*F)*W));
}

double H0f(r,B,X,m,p,A,H1,K,w)
double r,B,X,m,p,A,H1,K,w;
{
double H0=8*pi*r*r*r/sqrt(B)*X-(0.5*l*(l+1)*(m+4*pi*r*r*r*p)-w*w*r*r*r/A/B)*H1+(0.5*(l+2)*(l-1)*r-w*w*r*r*r/B-1/r*A*(m+4*pi*r*r*r*p)*(3*m-r+4*pi*r*r*r*p))*K;
H0=H0/(3*m+0.5*(l+2)*(l-1)*r+4*pi*r*r*r*p);
return(H0);
}

double Vf(r,w,e,p,B,A,Dp,W,H0,X)
double r,w,e,p,B,A,Dp,W,H0,X;
{
return(1/w/w/(e+p)*B*(1/sqrt(B)*X+1/r*Dp/sqrt(A)*W-0.5*(e+p)*H0));
}

double gf(e)
double e;
{
double le, p, gamma;
int i=j-1;
e=e/G*c*c;
le=log10(e);
for(;i>0;i--)
  {
  if(le>a[i-1])
    {
    p=b[i-1]*(le-a[i])/(a[i-1]-a[i])+b[i]*(le-a[i-1])/(a[i]-a[i-1]);
    p=pow(10,p);
    gamma=(e+p/c/c)/e*(b[i]-b[i-1])/(a[i]-a[i-1]);
    return(gamma);
    }
  }
return((pow(10,a[0])+pow(10,b[0])/c/c)/e*(b[i]-b[i-1])/(a[i]-a[i-1]));
}

double Ff(r,A,e,p,m,Dp)
double r,A,e,p,m,Dp;
{
double F;
F=8*pi/r/r*sqrt(A);
F=F*(e+3*p+r*Dp-(m/4/pi/r/r/r+p)*(4+A*(m/r-4*pi*r*r*e)));
return(F);
}

double getM(e0)
double e0;
{
	double dr=0.5,power=pow(10,b[0]),r=1,e,p,m,wpf=-1,p1,p2,p3,p4,m1,m2,m3,m4;
	e=e0;
	p=che(e);
	m=1.3333333*r*r*pi*e*r;
	for(;p>power;r=r+dr)
		{
			if(wpf>49995) return 0;
			wpf++;
			p1=fp(r,p,e,m);
			m1=fm(r,e);
			if((p+dr*p1/2)>power) e=ch(p+dr*p1/2); else break;
			p2=fp(r+dr/2,p+dr*p1/2,e,m+dr*m1/2);
			m2=fm(r+dr/2,e);
			if((p+dr*p2/2)>power) e=ch(p+dr*p2/2); else break;
			p3=fp(r+dr/2,p+dr*p2/2,e,m+dr*m2/2);
			m3=fm(r+dr/2,e);
			if((p+dr*p3)>power) e=ch(p+dr*p3); else break;
			p4=fp(r+dr,p+dr*p3,e,m+dr*m3);
			m4=fm(r+dr,e);
			p=p+dr*(p1+2*p2+2*p3+p4)/6;
			e=ch(p);
			m=m+dr*(m1+2*m2+2*m3+m4)/6;
		}
	return(m/M);
}


double getMmax()
{
	int n;
	double m0,m1,ma,mb,mc,md,E0,E1,Ea,Eb,Ec,Ed,Mmax,dE=1e12;

	E1=(E_core<4.2e18)?E_core:4.2e18;
	E0=E1-dE;
	m1=getM(E1);m0=getM(E0);
	if(m0-m1<=0)
	{
		Emax=E1;
		return m1;
	}
	for(n=0;n<7;n++)
	{
		E1=E1*0.72;
		E0=E1-dE;
		m1=getM(E1);m0=getM(E0);
		if(m0-m1<=0) break;
	}
	if(m0-m1>0)
	{
		Emax=E0;
		return m0;
	}
	Ea=E1;Eb=E1*1.3888888889;
	Ec=Ea+0.382*(Eb-Ea);
	Ed=Eb-Ec+Ea;
	mc=getM(Ec);md=getM(Ed);
	for(n=0;n<15;n++)
	{
		if(mc-md>0)
		{
			Eb=Ed;mb=md;
			Ed=Ec;md=mc;
			Ec=Ea+0.382*(Eb-Ea);
			mc=getM(Ec);
		}else{
			Ea=Ec;ma=mc;
			Ec=Ed;mc=md;
			Ed=Ea+0.618*(Eb-Ea);
			md=getM(Ed);
		}
	}
	Emax=(mc-md>0)?(0.5*(Ea+Ed)):(0.5*(Ec+Eb));
	return(getM(Emax));	
}

double getMmin()
{
	int n;
	double E0=7.2e17,m0;
	m0=getM(E0);
	for(n=0;n<7;n++)
	{
		if(m0<0.72) 
		{
			Emin=E0;
			return m0;
		}
		E0=E0*0.86;
		m0=getM(E0);
	}
	Emin=E0;
	return m0;
}

double M2rhoc(fM)
double fM;
{
	double Mmax,Mmin,E0,E1,E2,m0,m1,m2;
	int n;
	Mmax=getMmax();
	scr++;printf("\b\b\b\b\b\b 50%% %c",dysm[scr%6]);
	Mmin=getMmin();
	if(Mmax-fM<0) return Emax;
	if(Mmin-fM>0) return Emin;
	E0=Emin;E2=Emax;m0=Mmin;m2=Mmax;
	for(n=0;n<32;n++)
	{
		E1=0.5*(E0+E2);
		m1=getM(E1);
		if(m1-fM<0)
		{
			E0=E1;m0=m1;
		}else{
			E2=E1;m2=m1;
		}
	}
	return(E0+(E2-E0)*(fM-m0)/(m2-m0));
}

int fmode(e0)
double e0;
{
	double dr=0.5;
	double r,r0=1,R,RR,drx,rx;
	double ne,p,e,m,mR,A,B=1.0,BR,Bfactor,m1,m2,m3,m4,p1,p2,p3,p4,B1,B2,B3,B4,power,I=0,J,DDf,Df=0,f=1;
	double H1,H0,K,W,X,F,V,Dv,gamma,Dp,N;
	double DH11,DK1,DW1,DX1,H01,H02,K1,K2,x,X1,X2,Xp1,Xp2,W1,W2,V1,V2,V01,V02,V0;
	double w,wcheck,o[2],wi,aR,bR,gR,hR,kR,n,Y1,Y2,Z,DZ,DDZ,VZ,Ar1,Ar2,Ai1,Ai2,ar,ai,Ar[2],Ai[2],Br[2],Bi[2];
	int t,q,wpf,rpf,pfEnd,n4;
	scr++;printf("\b%c",dysm[scr%6]);
	power=pow(10,b[0]);
	r=r0;
	ne=e0;
	e=e0;
	p=che(e);
	m=1.3333333*r*r*pi*e*r;
	wpf=-1;
	for(;p>power;r=r+dr)
		{
			if(wpf>49995)
				{
					temp_rhoc=ne;
					temp_M=0;
					temp_R=0;
					temp_freq=0;
					temp_dpt=0;
					return 1;
				}
			wpf++;
			rhofile[wpf]=e*Gc2;   /*-- e,p,m in G=c=1 --*/
			pfile[wpf]=p*Gc4;
			Bfile[wpf]=B;
			mfile[wpf]=m*Gc2;
			A=1/(1-2*m*Gc2/r);
			p1=fp(r,p,e,m);
			m1=fm(r,e);
			if((p+dr*p1/2)>power) e=ch(p+dr*p1/2); else break;
			B1=Bf(r,p,m,B);
			p2=fp(r+dr/2,p+dr*p1/2,e,m+dr*m1/2);
			m2=fm(r+dr/2,e);
			if((p+dr*p2/2)>power) e=ch(p+dr*p2/2); else break;
			B2=Bf(r+dr/2,p+dr*p1/2,m+dr*m1/2,B+dr*B1/2);
			p3=fp(r+dr/2,p+dr*p2/2,e,m+dr*m2/2);
			m3=fm(r+dr/2,e);
			if((p+dr*p3)>power) e=ch(p+dr*p3); else break;
			B3=Bf(r+dr/2,p+dr*p2/2,m+dr*m2/2,B+dr*B2/2);
			p4=fp(r+dr,p+dr*p3,e,m+dr*m3);
			m4=fm(r+dr,e);
			B4=Bf(r+dr,p+dr*p3,m+dr*m3,B+dr*B3);
			J=-4*pi*(e+p/c/c)*Gc2*r*A;
			DDf=-(4/r*Df+J*Df+4/r*J*f);
			I=I-2.0/3*f*J/sqrt(A*B)*r*r*r*dr;
			p=p+dr*(p1+2*p2+2*p3+p4)/6;
			e=ch(p);
			m=m+dr*(m1+2*m2+2*m3+m4)/6;
			B=B+dr*(B1+2*B2+2*B3+B4)/6;
			f=f+Df*dr;
			Df=Df+DDf*dr;
		}
		pfEnd=wpf;
		R=r;
		mR=m*Gc2;
		BR=1-2*Gc2*m/r;
		Bfactor=BR/B;
		gamma=(b[1]-b[0])/(a[1]-a[0]);
		N=1/(gamma-1);
		RR=R-(N+1)*(p-dr*(p1+2*p2+2*p3+p4)/6)/(p1+2*p2+2*p3+p4)*6;
		wpf=-1;
		rpf=-1;
		scr++;printf("\b%c",dysm[scr%6]);
		for(n4=0;n4<pfEnd+1;n4++)
			{
				Bcor[n4]=Bfile[n4]*Bfactor;
			}
		I=I/sqrt(Bfactor);
		I=I/(f+2*I/r/r/r)/Gc2;
		I=m*sqrt(m/I)*Gc2;
		o[0]=(-0.0047+0.133*I+0.575*I*I)/mR-0.1e-5;
		o[1]=o[0]+0.2e-5;
		q=1;
		wcheck=0;
		for(t=0;;t++)
		{
			if(t>25)
			{
				temp_rhoc=ne;
				temp_M=0;
				temp_R=0;
				temp_freq=0;
				temp_dpt=0;
				return 2;
			}
			scr++;printf("\b%c",dysm[scr%6]);
			if(t==0) w=o[t];
			else w=o[q];
			e=rhofile[0];
			p=pfile[0];
			B=Bcor[0];
			W=1.0;
			K=(e+p);
			X=(e+p)*sqrt(B)*((4*pi/3*(e+3*p)-w*w/B/l)*W+0.5*K); 
			H1=(2*l*K+16*pi*(e+p)*W)/l/(l+1);
			rpf=-1;
			wpf=-1;
			r=r0;
			while(rpf<pfEnd)
			{
				rpf++;
				p=pfile[rpf];
				e=rhofile[rpf];
				B=Bcor[rpf];
				m=mfile[rpf];
				Dp=-(e+p)*(m+4*pi*r*r*r*p)/r/r/(1-2*m/r);
				Dv=-2*Dp/(e+p);
				A=1/(1-2*m/r);
				gamma=gf(e);
				H0=H0f(r,B,X,m,p,A,H1,K,w);
				V=Vf(r,w,e,p,B,A,Dp,W,H0,X);
				if(r==r0)V01=V;
				if(fabs((wcheck-w)/w)<Cri)
				{
					wpf++;
					Wfile[wpf]=sqrt(1-2*m/r)*W;
					Vfile[wpf]=V;
				}
				F=Ff(r,A,e,p,m,Dp);
				DH11=DH1(r,m,A,p,e,H1,H0,K,V);
				DK1=DK(r,H0,H1,Dv,K,e,p,A,W);
				DW1=DW(r,W,A,gamma,p,B,X,V,H0,K);
				DX1=DX(r,X,e,p,B,Dv,H0,w,H1,K,V,A,F,W);
				H1=H1+DH11*dr;
				K=K+DK1*dr;
				W=W+DW1*dr;
				X=X+DX1*dr;
				r=r+dr;
			}
			wpf=-1;
			rpf=-1;
			X1=X;Xp1=DX1;K1=K;H01=H0f(r,B,X,m,p,A,H1,K,w);W1=W;V1=V;
			scr++;printf("\b%c",dysm[scr%6]);
			p=pfile[0];
			e=rhofile[0];
			B=Bcor[0];
			W=1.0;K=-(e+p);
			X=(e+p)*sqrt(B)*((4*pi/3*(e+3*p)-w*w/B/l)*W+0.5*K);
			H1=(2*l*K+16*pi*(e+p)*W)/l/(l+1);
			r=r0;
			while(rpf<pfEnd)
			{
				rpf++;
				p=pfile[rpf];
				e=rhofile[rpf];
				B=Bcor[rpf];
				m=mfile[rpf];
				Dp=-(e+p)*(m+4*pi*r*r*r*p)/r/r/(1-2*m/r);
				Dv=-2*Dp/(e+p);
				A=1/(1-2*m/r);
				gamma=gf(e);
				H0=H0f(r,B,X,m,p,A,H1,K,w);
				V=Vf(r,w,e,p,B,A,Dp,W,H0,X);
				if(r==r0)V02=V;
				if((wcheck-w)/w>-Cri&&(wcheck-w)/w<Cri)
				{
					wpf++;
					Wfile[wpf]=sqrt(1-2*m/r)*W;
					Vfile[wpf]=V;
				}
				F=Ff(r,A,e,p,m,Dp);
				DH11=DH1(r,m,A,p,e,H1,H0,K,V);
				DK1=DK(r,H0,H1,Dv,K,e,p,A,W);
				DW1=DW(r,W,A,gamma,p,B,X,V,H0,K);
				DX1=DX(r,X,e,p,B,Dv,H0,w,H1,K,V,A,F,W);
				H1=H1+DH11*dr;
				K=K+DK1*dr;
				W=W+DW1*dr;
				X=X+DX1*dr;
				r=r+dr;
			}
			wpf=-1;
			rpf=-1;
			X2=X;Xp2=DX1;K2=K;H02=H0f(r,B,X,m,p,A,H1,K,w);W2=W;V2=V;
			x=-(X1-(RR-R)/(N+1)*Xp1)/(X2-(RR-R)/(N+1)*Xp2);
			H0=H01+x*H02; K=K1+x*K2; W=W1+x*W2; V=V1+x*V2;V0=V01+x*V02;
			if(fabs((wcheck-w)/w)<Cri)
			{
				r=r0;
				while(rpf<pfEnd)
				{
					rpf++;
					W1=Wfile1[rpf];
					W2=Wfile2[rpf];
					W=W1+x*W2;
					wpf++;
					Wfile[wpf]=W/(1+x);
					V1=Vfile1[rpf];
					V2=Vfile2[rpf];
					V=V1+x*V2;
					Vfile[wpf]=V/V0;
					r=r+dr;
				}
				break;
			}
			wcheck=w;
			n=0.5*(l-1)*(l+2);
			aR=-(n*R+3*mR)/(w*w*R*R-(n+1)*mR/R);
			bR=(n*R*(R-2*mR)-w*w*R*R*R*R+mR*(R-3*mR));
			bR=bR/(R-2*mR)/(w*w*R*R-(n+1)*mR/R);
			gR=n*(n+1)*R*R+3*n*mR*R+6*mR*mR;
			gR=gR/R/R/(n*R+3*mR);
			hR=-n*R*R+3*n*mR*R+3*mR*mR;
			hR=hR/(R-2*mR)/(n*R+3*mR);
			kR=-R*R/(R-2*mR);
			Y1=K;
			Y2=aR*H0+bR*K;
			Z=(kR*Y1-Y2)/(kR*gR-hR);
			DZ=(gR*Y2-hR*Y1)/(gR*kR-hR);
			if(w<5e-9)
				{
					temp_rhoc=ne;
					temp_M=0;
					temp_R=0;
					temp_freq=0;
					temp_dpt=0;
					return 3;
				}
			for(r=R;r<25.0/w;r=r+dr)
			{
				drx=dr/(1-2*mR/r);
				VZ=(1-2*mR/r)/r/r/r/(n*r+3*mR)/(n*r+3*mR);
				VZ=VZ*(2.0*n*n*(n+1)*r*r*r+6.0*n*n*mR*r*r+18.0*n*mR*mR*r+18*mR*mR*mR);
				DDZ=(VZ-w*w)*Z;
				Z=Z+DZ*drx;
				DZ=DZ+DDZ*drx;
			}
			r=r-dr;
			rx=r+2*mR*log(r/2/mR-1);
			Ar1=2*cos(w*rx)-2*(n+1)/w/r*sin(w*rx)+1/w/w/r/r*(1.5*mR*w*(1+2/n)*sin(w*rx)-n*(n+1)*cos(w*rx));
			Ai1=2*sin(w*rx)+2*(n+1)/w/r*cos(w*rx)-1/w/w/r/r*(1.5*mR*w*(1+2/n)*cos(w*rx)+n*(n+1)*sin(w*rx));
			Ar2=-2*w*sin(w*rx)-2*(n+1)*cos(w*rx)/r+1/w/r/r*(1.5*mR*w*(1+2/n)*cos(w*rx)+n*(n+1)*sin(w*rx))+(1-2*mR/r)*2*(n+1)/w/r/r*sin(w*rx);
			Ai2=2*w*cos(w*rx)-2*(n+1)*sin(w*rx)/r+1/w/r/r*(1.5*mR*w*(1+2/n)*sin(w*rx)-n*(n+1)*cos(w*rx))-(1-2*mR/r)*2*(n+1)/w/r/r*cos(w*rx);
			ar=(Ai2*Z-Ai1*DZ)/(Ar1*Ai2-Ar2*Ai1);
			ai=-(Ar1*DZ-Ar2*Z)/(Ar1*Ai2-Ar2*Ai1);
			if(t==0)
			{
				Ar[t]=ar;
				Ai[t]=ai;
			}else{
				Ar[q]=ar;
				Ai[q]=ai;
				Br[0]=(o[0]*Ar[1]-o[1]*Ar[0])/(o[0]-o[1]);
				Br[1]=(Ar[0]-Ar[1])/(o[0]-o[1]);
				Bi[0]=(o[0]*Ai[1]-o[1]*Ai[0])/(o[0]-o[1]);
				Bi[1]=(Ai[0]-Ai[1])/(o[0]-o[1]);
				w=-(Br[0]*Br[1]+Bi[0]*Bi[1])/(Br[1]*Br[1]+Bi[1]*Bi[1]);
				if (w<=o[0])
				{
					o[1]=o[0];o[0]=w;Ar[1]=Ar[0];Ai[1]=Ai[0];q=0;
				}else if(w>=o[1]){
					o[0]=o[1];o[1]=w;Ar[0]=Ar[1];Ai[0]=Ai[1];q=1;
				}else if((o[1]-w)>(w-o[0])){
					o[1]=w;q=1;
				}else{
					o[0]=w;q=0;
				}
			}
		}
		wi=(Br[0]*Bi[1]-Bi[0]*Br[1])/(Br[1]*Br[1]+Bi[1]*Bi[1]);
		temp_rhoc=ne;
		temp_M=mR;
		temp_M=temp_M/(M*Gc2);
		temp_R=RR/1000;
		temp_freq=w*c/(2000*pi);
		temp_dpt=1/(wi*c);
		scr++;printf("\b%c",dysm[scr%6]);
	return 0;
}

int main()
{
	FILE *temp,*Xout;
	int pf,n,len1,len2,neos,mode;
	double srhoc[40],smass[40],sradius[40],sfreq[40],sdptime[40];
	char pathName[50];
	time_t timep;
	struct tm *sysTime;
	time(&timep);
	sysTime=gmtime(&timep);
	system("title Advanced codes F-mode_v37     for Windows");
	temp=fopen("temp_dir.dat","w");
	fclose(temp);
	remove("temp_dir.dat");
	system("dir EoS_lib /B /ON /A-D-R-H >> temp_dir.dat");
	if((temp=fopen("temp_dir.dat","r"))==NULL){printf("No Equation of State files detected, press any key to leave.");fclose(temp);system("pause");return 0;
	}else{printf("Equation of State files detected:\n\n");}
	n=0;
	while(fscanf(temp,"%s",list[n].Name)>0){printf("%d\t%s\n",n+1,list[n].Name);n++;}
	len1=n;
	fclose(temp);
	remove("temp_dir.dat");
	printf("\nCaution: please ensure all *.txt files above are available EoS.\n\nPress any key to start (Ctrl+C to quit).\n");
	mfile=(double*)malloc(50000*sizeof(double));
	pfile=(double*)malloc(50000*sizeof(double));
	rhofile=(double*)malloc(50000*sizeof(double));
	Bfile=(double*)malloc(50000*sizeof(double));
	Bcor=(double*)malloc(50000*sizeof(double));
	Wfile=(double*)malloc(50000*sizeof(double));
	Wfile1=(double*)malloc(50000*sizeof(double));
	Wfile2=(double*)malloc(50000*sizeof(double));
	Vfile=(double*)malloc(50000*sizeof(double));
	Vfile1=(double*)malloc(50000*sizeof(double));
	Vfile2=(double*)malloc(50000*sizeof(double));
	system("pause");
	printf("\nPlease select a work mode:\n\t[0]\tCompute within a density range\n\t[1]\tSearch for a given mass\n\t[2]\tSearch for a given frequency\n\nEnter the mode number: ");
	fflush(stdin);
	scanf("%d",&mode);
	if(mode==0)
	{
		double strhoc,edrhoc,deltrho;
		printf("\nTips: 3e17<rhoc<3e18\n\nStart rhoc(kg.m-3)=");
		fflush(stdin);
	    scanf("%lf",&strhoc);
		while((strhoc<3e17)||(strhoc>3e18))
		{
			printf("\nWrong input, Start at rhoc(kg.m-3)=");
		    fflush(stdin);
	        scanf("%lf",&strhoc);
		}
		printf("\nEnd at rhoc(kg.m-3)=");
		fflush(stdin);
	    scanf("%lf",&edrhoc);
		while((edrhoc<3e17)||(edrhoc>3e18))
		{
			printf("\nWrong input, End at rhoc(kg.m-3)=");
		    fflush(stdin);
	        scanf("%lf",&edrhoc);
		}
		outf=fopen("Results_mod_0.txt","a");
		Xout=fopen("Results_mod_0.xls","a");
		fprintf(outf,"\nTime: %d-%d-%d-%dh-%dmin-%ds\n",1900+sysTime->tm_year,1+sysTime->tm_mon,sysTime->tm_mday,8+sysTime->tm_hour,sysTime->tm_min,sysTime->tm_sec);
		fprintf(outf,"====================================\n");
		fprintf(outf,"%-25s\t%-15s\t%-10s\t%-10s\t%-15s\t%-18s\n","EoS","rhoc(kg.m-3)","M(Msun)","R(km)","frequency(kHz)","dampingTime(s)");
		fprintf(Xout,"\nTime: %d-%d-%d-%dh-%dmin-%ds\n",1900+sysTime->tm_year,1+sysTime->tm_mon,sysTime->tm_mday,8+sysTime->tm_hour,sysTime->tm_min,sysTime->tm_sec);
		fprintf(Xout,"====================================\n");
		fprintf(Xout,"%-25s\t%-15s\t%-10s\t%-10s\t%-15s\t%-18s\n","EoS","rhoc(kg.m-3)","M(Msun)","R(km)","frequency(kHz)","dampingTime(s)");
		for(neos=1;neos<=len1;neos++)
		{
			j=0;scr=0;
			sprintf(pathName,".\\EoS_lib\\%s",list[neos-1].Name);
			if((inf=fopen(pathName,"r"))==NULL)
				{printf("\nCannot load %s!\n",list[neos-1].Name);continue;}
			while(fscanf(inf,"%lf",&a[j])==1){fscanf(inf,"%lf%*[^\n]",&b[j]);j++;}
			fclose(inf);
			len2=strlen(list[neos-1].Name);
			for(pf=0;pf<=len2-2;pf++)
			{
				if(list[neos-1].Name[pf]=='.'){break;}
				list[neos-1].tName[pf]=list[neos-1].Name[pf];
			}
			if(strhoc==edrhoc)
			{
				printf("\n%s is under computed...      %c",list[neos-1].tName,dysm[scr%6]);
				fmode(strhoc);
				printf("\b \n%-8s\t%-8s\t%-15s\t%-18s\n","M(Msun)","R(km)","frequency(kHz)","dampingTime(s)");
				printf("%-8.2lf\t%-8.2lf\t%-15lf\t%-18lf\n",temp_M,temp_R,temp_freq,temp_dpt);
				fprintf(outf,"\n%-25s\t%-15e\t%-10lf\t%-10lf\t%-15lf\t%-18lf\n",list[neos-1].tName,strhoc,temp_M,temp_R,temp_freq,temp_dpt);
				fprintf(Xout,"\n%-25s\t%-15e\t%-10lf\t%-10lf\t%-15lf\t%-18lf\n",list[neos-1].tName,strhoc,temp_M,temp_R,temp_freq,temp_dpt);
				continue;
			}
			if(edrhoc>strhoc)
			{
				edrhoc=(edrhoc<pow(10,a[j-1]))?edrhoc:pow(10,a[j-1]);
			}else{
				strhoc=(strhoc<pow(10,a[j-1]))?strhoc:pow(10,a[j-1]);
			}
			deltrho=0.03125*(edrhoc-strhoc);
			printf("\n%s is under computed...      %c",list[neos-1].tName,dysm[scr%6]);
			for(pf=0;pf<33;pf++)
			{
				srhoc[pf]=strhoc+pf*deltrho;
				fmode(srhoc[pf]);
				smass[pf]=temp_M;sradius[pf]=temp_R;sfreq[pf]=temp_freq;sdptime[pf]=temp_dpt;
				scr++;printf("\b\b\b\b\b\b%3d%% %c",pf*3,dysm[scr%6]);
				fprintf(outf,"\n%-25s\t%-15e\t%-10lf\t%-10lf\t%-15lf\t%-18lf\n",list[neos-1].tName,srhoc[pf],smass[pf],sradius[pf],sfreq[pf],sdptime[pf]);
				fprintf(Xout,"\n%-25s\t%-15e\t%-10lf\t%-10lf\t%-15lf\t%-18lf\n",list[neos-1].tName,srhoc[pf],smass[pf],sradius[pf],sfreq[pf],sdptime[pf]);
			}
			scr++;printf("\b\b\b\b\b\b100%% %c",dysm[scr%6]);
			printf("\b ");
		}
		fclose(outf);
		fclose(Xout);
	}else if(mode==1)
	{
		double objM,objrho;
		printf("\nTips: 1.0<M<2.74\n\nM(Msun)=");
		fflush(stdin);
	    scanf("%lf",&objM);
		while((objM<0.5)||(objM>3.0))
		{
			printf("\nWrong input, M(Msun)=");
		    fflush(stdin);
	        scanf("%lf",&objM);
		}
		outf=fopen("Results_mod_1.txt","a");
		Xout=fopen("Results_mod_1.xls","a");
		fprintf(outf,"\nTime: %d-%d-%d-%dh-%dmin-%ds\n",1900+sysTime->tm_year,1+sysTime->tm_mon,sysTime->tm_mday,8+sysTime->tm_hour,sysTime->tm_min,sysTime->tm_sec);
		fprintf(outf,"====================================\n");
		fprintf(outf,"%-25s\t%-15s\t%-10s\t%-10s\t%-15s\t%-18s\n","EoS","rhoc(kg.m-3)","M(Msun)","R(km)","frequency(kHz)","dampingTime(s)");
		fprintf(Xout,"\nTime: %d-%d-%d-%dh-%dmin-%ds\n",1900+sysTime->tm_year,1+sysTime->tm_mon,sysTime->tm_mday,8+sysTime->tm_hour,sysTime->tm_min,sysTime->tm_sec);
		fprintf(Xout,"====================================\n");
		fprintf(Xout,"%-25s\t%-15s\t%-10s\t%-10s\t%-15s\t%-18s\n","EoS","rhoc(kg.m-3)","M(Msun)","R(km)","frequency(kHz)","dampingTime(s)");
		for(neos=1;neos<=len1;neos++)
		{
			j=0;scr=0;
			sprintf(pathName,".\\EoS_lib\\%s",list[neos-1].Name);
			if((inf=fopen(pathName,"r"))==NULL)
				{printf("\nCannot load %s!\n",list[neos-1].Name);continue;}
			while(fscanf(inf,"%lf",&a[j])==1){fscanf(inf,"%lf%*[^\n]",&b[j]);j++;}
			E_core=pow(10,a[j-1]);
			fclose(inf);
			len2=strlen(list[neos-1].Name);
			for(pf=0;pf<=len2-2;pf++)
			{
				if(list[neos-1].Name[pf]=='.'){break;}
				list[neos-1].tName[pf]=list[neos-1].Name[pf];
			}
			scr=scr%6;printf("\n%s is under computed...      %c",list[neos-1].tName,dysm[scr%6]);
			scr++;printf("\b\b\b\b\b\b  0%% %c",dysm[scr%6]);
			objrho=M2rhoc(objM);
			fmode(objrho);
			scr++;printf("\b\b\b\b\b\b100%% %c",dysm[scr%6]);
			printf("\b f_%.2lf=%.3lf(kHz), dT=%.3lf(s)\n",temp_M,temp_freq,temp_dpt);
			n=pf;
			fprintf(outf,"\n%-25s\t%-15e\t%-10lf\t%-10lf\t%-15lf\t%-18lf\n",list[neos-1].tName,objrho,temp_M,temp_R,temp_freq,temp_dpt);
			fprintf(Xout,"\n%-25s\t%-15e\t%-10lf\t%-10lf\t%-15lf\t%-18lf\n",list[neos-1].tName,objrho,temp_M,temp_R,temp_freq,temp_dpt);
		}
		fclose(outf);
		fclose(Xout);
    }else if(mode==2)
	{
		double objF,lefrhoc,midrhoc,rigrhoc,lefF,midF,rigF,Mmax,Mmin,lefT,rigT,lefR,rigR;
		printf("\nTips: 0.5<freq<3.6\n\nfreq(kHz)=");
		fflush(stdin);
	    scanf("%lf",&objF);
		while((objF<0.5)||(objF>3.6))
		{
			printf("\nWrong input, freq(kHz)=");
		    fflush(stdin);
	        scanf("%lf",&objF);
		}
		outf=fopen("Results_mod_2.txt","a");
		Xout=fopen("Results_mod_2.xls","a");
		fprintf(outf,"\nTime: %d-%d-%d-%dh-%dmin-%ds\n",1900+sysTime->tm_year,1+sysTime->tm_mon,sysTime->tm_mday,8+sysTime->tm_hour,sysTime->tm_min,sysTime->tm_sec);
		fprintf(outf,"====================================\n");
		fprintf(outf,"%-25s\t%-15s\t%-10s\t%-10s\t%-15s\t%-18s\n","EoS","rhoc(kg.m-3)","M(Msun)","R(km)","frequency(kHz)","dampingTime(s)");
		fprintf(Xout,"\nTime: %d-%d-%d-%dh-%dmin-%ds\n",1900+sysTime->tm_year,1+sysTime->tm_mon,sysTime->tm_mday,8+sysTime->tm_hour,sysTime->tm_min,sysTime->tm_sec);
		fprintf(Xout,"====================================\n");
		fprintf(Xout,"%-25s\t%-15s\t%-10s\t%-10s\t%-15s\t%-18s\n","EoS","rhoc(kg.m-3)","M(Msun)","R(km)","frequency(kHz)","dampingTime(s)");
		for(neos=1;neos<=len1;neos++)
		{
			j=0;scr=0;
			sprintf(pathName,".\\EoS_lib\\%s",list[neos-1].Name);
			if((inf=fopen(pathName,"r"))==NULL)
				{printf("\nCannot load %s!\n",list[neos-1].Name);continue;}
			while(fscanf(inf,"%lf",&a[j])==1){fscanf(inf,"%lf%*[^\n]",&b[j]);j++;}
			fclose(inf);
			E_core=pow(10,a[j-1]);
			len2=strlen(list[neos-1].Name);
			for(pf=0;pf<=len2-2;pf++)
			{
				if(list[neos-1].Name[pf]=='.'){break;}
				list[neos-1].tName[pf]=list[neos-1].Name[pf];
			}
			scr=scr%6;printf("\n%s is under computed...      %c",list[neos-1].tName,dysm[scr%6]);
			Mmin=getMmin();
			lefrhoc=Emin;
			Mmax=getMmax();
			rigrhoc=Emax;
			fmode(lefrhoc);
			lefF=temp_freq;
			lefT=temp_dpt;
			lefR=temp_R;
			fmode(rigrhoc);
			rigF=temp_freq;
			rigT=temp_dpt;
			rigR=temp_R;
			scr++;printf("\b\b\b\b\b\b  0%% %c",dysm[scr%6]);
			if((lefF<objF)&&(objF<rigF))
		    {
			for(pf=0;pf<32;pf++)
			{
				scr++;printf("\b\b\b\b\b\b%3d%% %c",pf*3,dysm[scr%6]);
				midrhoc=0.5*(lefrhoc+rigrhoc);
				fmode(midrhoc);
				midF=temp_freq;
				if(midF<objF)
					{
						lefrhoc=midrhoc;
						lefF=midF;
					}else{
						rigrhoc=midrhoc;
						rigF=midF;
					}
				
			}
			midrhoc=0.5*(lefrhoc+rigrhoc);
			fmode(midrhoc);
			scr++;printf("\b\b\b\b\b\b100%% %c",dysm[scr%6]);
			printf("\b f_%.2lf=%.3lf(kHz), dT=%.3lf(s)\n",temp_M,temp_freq,temp_dpt);
			fprintf(outf,"\n%-25s\t%-15e\t%-10lf\t%-10lf\t%-15lf\t%-18lf\n",list[neos-1].tName,midrhoc,temp_M,temp_R,temp_freq,temp_dpt);
			fprintf(Xout,"\n%-25s\t%-15e\t%-10lf\t%-10lf\t%-15lf\t%-18lf\n",list[neos-1].tName,midrhoc,temp_M,temp_R,temp_freq,temp_dpt);
			}else if(objF>=rigF){
				scr++;printf("\b\b\b\b\b\b100%% %c",dysm[scr%6]);
				printf("\b f_%.2lf=%.3lf(kHz), dT=%.3lf(s)\n",Mmax,rigF,rigT);
				fprintf(outf,"\n%-25s\t%-15e\t%-10lf\t%-10lf\t%-15lf\t%-18lf\n",list[neos-1].tName,rigrhoc,Mmax,rigR,rigF,rigT);
				fprintf(Xout,"\n%-25s\t%-15e\t%-10lf\t%-10lf\t%-15lf\t%-18lf\n",list[neos-1].tName,rigrhoc,Mmax,rigR,rigF,rigT);
			}else{
				scr++;printf("\b\b\b\b\b\b100%% %c",dysm[scr%6]);
				printf("\b f_%.2lf=%.3lf(kHz), dT=%.3lf(s)\n",Mmin,lefF,lefT);
				fprintf(outf,"\n%-25s\t%-15e\t%-10lf\t%-10lf\t%-15lf\t%-18lf\n",list[neos-1].tName,lefrhoc,Mmin,lefR,lefF,lefT);
				fprintf(Xout,"\n%-25s\t%-15e\t%-10lf\t%-10lf\t%-15lf\t%-18lf\n",list[neos-1].tName,lefrhoc,Mmin,lefR,lefF,lefT);
			}
		}
		fclose(outf);
		fclose(Xout);
	}
	free(mfile);
	free(pfile);
	free(rhofile);
	free(Bfile);
	free(Bcor);
	free(Wfile);
	free(Wfile1);
	free(Wfile2);
	free(Vfile);
	free(Vfile1);
	free(Vfile2);
	system("pause");
	return 0;
}