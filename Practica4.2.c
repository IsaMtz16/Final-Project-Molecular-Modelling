#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 125
#define Nt 100000

extern float Random(void);
#define NormRANu (2.3283063671E-10)
extern void ini_ran(int SEMILLA);
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran, ig1, ig2, ig3;

#define PI 3.14159265

typedef struct{
    double x, y, z;
} Vector;

Vector particula[N]={0};
Vector velocidad[N]={0};


float Random(void)    // generar numero aleatorio
{
    float r;
    ig1=ind_ran-24;
    ig2=ind_ran-55;
    ig3=ind_ran-61;
    irr[ind_ran]=irr[ig1]+irr[ig2];
    ir1=(irr[ind_ran]^irr[ig3]);
    ind_ran++;
    r=ir1*NormRANu;
    return r;
}

void ini_ran(int SEMILLA)
{
    int INI,FACTOR,SUM,i;
    srand(SEMILLA);
    INI=SEMILLA;
    FACTOR=67397;
    SUM=7364893;
    for(i=0; i<256; i++)
    {
        INI=(INI*FACTOR+SUM);
        irr[i]=INI;
    }
    ind_ran=ig1=ig2=ig3=0;
}

/*void leer_fichero(Vector *particula){
    int  i;
    FILE*f;
    f=fopen("Conf_final.txt","r");
    for(i=0;i<N;i++){
        fscanf(f,"%lf %lf  %lf",&particula[i].x,&particula[i].y,&particula[i].z);
        //printf("%lf     %lf     %lf\n",particula[i].x,particula[i].y,particula[i].z);
    }
    fclose(f);

}*/

void inicializar_sc(Vector *particula,double *L,double rho){

    double M;
    double a;
    double exponente;

    exponente=(double)1/(double)3;

    M=pow(N,exponente);
    *L=pow(N/rho,exponente);
    a=*L/M;
    //printf("%lf",a);

    int i,j,k,n;

    //FICHERO PARA ESCRIBIR LA CONFIGURACIÓN INICIAL
    /*FILE*f;
    f=fopen("Conf_inicial_sc.xyz","w");
    fprintf(f,"%d\n\n",N);*/

    n=0;
    for(k=0;(k<M)&&(n<N);k++){ //Bucle para z
        for(j=0;(j<M)&&(n<N);j++){ //Bucle para y
            for(i=0;(i<M)&&(n<N);i++){ //Buclue para x
                    particula[n].x=(double)i*a;
                    particula[n].y=(double)j*a;
                    particula[n].z=(double)k*a;

                    //fprintf(f,"S    %lf     %lf     %lf\n",particula[n].x,particula[n].y,particula[n].z);
                    //printf("%lf\n",a);
                    n++;
            }
        }

    }

}

void inicializar_velocidad(Vector *velocidad,double T){
    int i;
    double a,b,c;
    for(i=0;i<N;i++){
        a=Random();
        b=Random();
        c=Random();
        if(a<0.5){
            velocidad[i].x=sqrt(T);
        }
        else{
            velocidad[i].x=-sqrt(T);
        }
        if(b<0.5){
            velocidad[i].y=sqrt(T);
        }
        else{
            velocidad[i].y=-sqrt(T);
        }
        if(c<0.5){
            velocidad[i].z=sqrt(T);
        }
        else{
            velocidad[i].z=-sqrt(T);
        }
    }

}

void Box_Muller (double *g1, double *g2,double sigma)
{
    double d1, d2;
    d1=sigma*sqrt(-2.0*log(Random()));
    d2=2.0*PI*Random();
    *g1=-d1*cos(d2);
    *g2=-d1*sin(d2);
}


void Anderson_thermostat(Vector *velocidad,double T,double nu){
    double sigma;
    sigma=sqrt(T);
    double a,g1,g2;
    int i;
    for(i=0;i<N;i++){
        a=Random();
        //printf("%lf     %lf     %lf\n",velocidad[i].x,velocidad[i].y,velocidad[i].z);
        if(a<nu){
            Box_Muller(&g1,&g2,sigma);
            velocidad[i].x=g1;
            velocidad[i].y=g2;
            Box_Muller(&g1,&g2,sigma);
            velocidad[i].z=g1;

        }
    }

}

void pbc(double *x,double L){
    if(*x>(L/2)){
        *x=*x-L;
    }
    if(*x<-(L/2)){
        *x=*x+L;
    }

}

void vector_0(Vector *fuerza){
    int i;
    for(i=0;i<N;i++){
        fuerza[i].x=0;
        fuerza[i].y=0;
        fuerza[i].z=0;

    }
}

void fuerza_energia(double L,Vector *particula,double cutoff,Vector *fuerza,double *Epot){
    *Epot=0;
    vector_0(fuerza);

    double dx, dy, dz,d;
    int i,j;
    for(i=0;i<N;i++){

        for(j=i+1;j<N;j++){
                dx=particula[i].x-particula[j].x;
                dy=particula[i].y-particula[j].y;
                dz=particula[i].z-particula[j].z;
                pbc(&dx,L);
                pbc(&dy,L);
                pbc(&dz,L);
                d=sqrt(dx*dx+dy*dy+dz*dz);

                //printf("%lf\n",d);

                if(d<cutoff){
                    fuerza[i].x=fuerza[i].x+(48/(pow(d,14))-24/(pow(d,8)))*dx;
                    fuerza[i].y=fuerza[i].y+(48/(pow(d,14))-24/(pow(d,8)))*dy;
                    fuerza[i].z=fuerza[i].z+(48/(pow(d,14))-24/(pow(d,8)))*dz;
                    //printf("%lf",48/(pow(d,14)));
                    //printf("%lf     %lf     %lf\n",fuerza[i].x,fuerza[i].y,fuerza[i].z);

                    fuerza[j].x=fuerza[j].x-(48/(pow(d,14))-24/(pow(d,8)))*dx;
                    fuerza[j].y=fuerza[j].y-(48/(pow(d,14))-24/(pow(d,8)))*dy;
                    fuerza[j].z=fuerza[j].z-(48/(pow(d,14))-24/(pow(d,8)))*dz;

                    *Epot=*Epot+(double)4*(1/(pow(d,12))-1/(pow(d,6)))-(double)4*(1/(pow(cutoff,12))-1/(pow(cutoff,6)));
                }


        }
    }

}
void Presion(double L,Vector Particula[],double *P,double rho,double T,double cutoff){
    int i,j;
    double dx, dy, dz,d;
    double aux=0;
    for(i=0;i<N;i++){
            for(j=i+1;j<N;j++){
                dx=particula[i].x-particula[j].x;
                dy=particula[i].y-particula[j].y;
                dz=particula[i].z-particula[j].z;
                pbc(&dx,L);
                pbc(&dy,L);
                pbc(&dz,L);
                d=sqrt(dx*dx+dy*dy+dz*dz);
                if(d<cutoff){
                    aux=aux+(48/(pow(d,14))-24/(pow(d,8)))*dx*dx+(48/(pow(d,14))-24/(pow(d,8)))*dy*dy+(48/(pow(d,14))-24/(pow(d,8)))*dz*dz;
                }

            }

    }
    aux=aux/((double)3*L*L*L);
    *P=T*rho+aux;

}

void Energia_cinetica(Vector velocidad[],double *Ecin){
    *Ecin=0;
    int i;
    for(i=0;i<N;i++){
        *Ecin=*Ecin+velocidad[i].x*velocidad[i].x+velocidad[i].y*velocidad[i].y+velocidad[i].z*velocidad[i].z;
    }
    *Ecin=*Ecin/2;
}

void Momento(Vector velocidad[],double *p){
    int i;
    *p=0;
    for(i=0;i<N;i++){
        *p=*p+sqrt(velocidad[i].x*velocidad[i].x+velocidad[i].y*velocidad[i].y+velocidad[i].z*velocidad[i].z);
    }
}

void Temperatura(double Ecin,double *T_inst){
        *T_inst=2*Ecin/(3*(double)N-3);
}

void modulo(Vector velocidad[],int i,double *v){
    *v=sqrt(velocidad[i].x*velocidad[i].x+velocidad[i].y*velocidad[i].y+velocidad[i].z*velocidad[i].z);
}

double reescalar(double E,double epsilon){
    return E*epsilon;
}

double reescalar_P(double P,double epsilon, double sigma){
    double epsilon_=epsilon*1000;
    double sigma_=sigma/100;
    return P*epsilon_/(sigma_*sigma_*sigma_*6.022*pow(10,23));
}



void velocity_Verlet(double h,double L,Vector *particula,double cutoff,double T,double nu,double *media_Epot,double *media_Ecin,double rho,double *media_P){
    int i,j,k,m;
    double t=0;
    double T_inst;
    double E, Epot, Ecin;
    double P;

    double p;//momento
    E=0;
    double v;//modulo de la velocidad

    Vector fuerza[N]={0};
    vector_0(velocidad);

    *media_Epot=0; *media_Ecin=0; *media_P=0;
    for(i=0;i<Nt;i++){

        t=t+h;
        fuerza_energia(L,particula,cutoff,fuerza,&Epot);
        Energia_cinetica(velocidad,&Ecin);


        E=Epot+Ecin;
        //Momento(velocidad,&p);


        for(j=0;j<N;j++){

            particula[j].x=particula[j].x+velocidad[j].x*h+0.5*fuerza[j].x*h*h;
            particula[j].y=particula[j].y+velocidad[j].y*h+0.5*fuerza[j].y*h*h;
            particula[j].z=particula[j].z+velocidad[j].z*h+0.5*fuerza[j].z*h*h;

             //PERIODIC BOUNDARY CONDITIONS
            pbc(&particula[j].x,L);
            pbc(&particula[j].y,L);
            pbc(&particula[j].z,L);

            velocidad[j].x=velocidad[j].x+0.5*fuerza[j].x*h;
            velocidad[j].y=velocidad[j].y+0.5*fuerza[j].y*h;
            velocidad[j].z=velocidad[j].z+0.5*fuerza[j].z*h;


        }
        fuerza_energia(L,particula,cutoff,fuerza,&Epot);
        for(k=0;k<N;k++){

            velocidad[k].x=velocidad[k].x+0.5*fuerza[k].x*h;
            velocidad[k].y=velocidad[k].y+0.5*fuerza[k].y*h;
            velocidad[k].z=velocidad[k].z+0.5*fuerza[k].z*h;


            Anderson_thermostat(velocidad,T,nu);
        }
        Energia_cinetica(velocidad,&Ecin);

        E=Epot+Ecin;
        *media_Epot=*media_Epot+Epot;
        *media_Ecin=*media_Ecin+Ecin;
        //Temperatura(Ecin,&T_inst);
        //Momento(velocidad,&p);

        Presion(L,particula,&P,rho,T,cutoff);
        *media_P=*media_P+P;
    }
    *media_Epot=*media_Epot/(double)Nt;
    *media_Ecin=*media_Ecin/(double)Nt;
    *media_P=*media_P/(double)Nt;

}



void fichero(double valores_rho[],double E[],double Epot[],double Ecin[]){
    FILE*g;
    g=fopen("Energia_densidad.txt","w");

    int i;
    for(i=0;i<6;i++){
        fprintf(g,"%lf      %lf     %lf     %lf\n",valores_rho[i],E[i],Epot[i],Ecin[i]);
    }
    fclose(g);
}

void fichero_P(double valores_rho[],double P[]){
    FILE*h;
    h=fopen("Presion.txt","w");

    int i;
    for(i=0;i<6;i++){
        fprintf(h,"%lf      %lf\n",valores_rho[i],P[i]);
    }
    fclose(h);
}

int main(){
    ini_ran(time(NULL));
    double cutoff=2.5;
    double T;
    double nu=0.1;
    double h=0.001;
    double L;

    double rho;

    double epsilon=0.998;//kJ/mol
    double sigma=3.4/(pow(10,8));//cm
    double m=40;//g/mol

    T=1.2;//T'

    double valores_rho[6] = {0.05,0.1,0.2,0.4,0.6,0.8};//rho=rho'*m/sigma^3
    double media_Epot, media_Ecin,media_P;
    double E[6]={0};
    double Epot[6]={0};
    double Ecin[6]={0};
    double P[6]={0};

    int i;
    for(i=0;i<6;i++){
        rho=valores_rho[i];
        printf("%lf\n",rho);

        inicializar_sc(particula,&L,rho);

        //inicializar_velocidad(velocidad,T);
        velocity_Verlet(h,L,particula,cutoff,T,nu,&media_Epot,&media_Ecin,rho,&media_P);

        Epot[i]=reescalar(media_Epot,epsilon);
        Ecin[i]=reescalar(media_Ecin,epsilon);
        E[i]=Epot[i]+Ecin[i];
        printf("%lf\n",media_P);
        P[i]=reescalar_P(media_P,epsilon,sigma);
        printf("%lf\n",P[i]);

        valores_rho[i]=valores_rho[i]*m/(pow(sigma,3));
    }
    //fichero(valores_rho,E,Epot,Ecin);
    fichero_P(valores_rho,P);

return 0;
}
