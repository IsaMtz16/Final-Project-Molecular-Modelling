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

void inicializar_sc(Vector *particula,double *L){
    double rho=0.7;

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
    FILE*f;
    f=fopen("Conf_inicial_sc.xyz","w");
    fprintf(f,"%d\n\n",N);

    n=0;
    for(k=0;(k<M)&&(n<N);k++){ //Bucle para z
        for(j=0;(j<M)&&(n<N);j++){ //Bucle para y
            for(i=0;(i<M)&&(n<N);i++){ //Buclue para x
                    particula[n].x=(double)i*a;
                    particula[n].y=(double)j*a;
                    particula[n].z=(double)k*a;

                    fprintf(f,"S    %lf     %lf     %lf\n",particula[n].x,particula[n].y,particula[n].z);
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

void Euler(Vector *particula, double L,double cutoff,double h){

    double E, Epot, Ecin;
    double p;
    Vector fuerza[N]={0};

    int i,t;
    double tiempo=0;
    FILE*f2;
    f2=fopen("Energia_momento_t_Euler_h_0.0001.txt","w");

    //ABRIMOS UN FICHERO PARA ESCRIBIR LAS ENERGÍAS
    //FILE*f3;
    //f3=fopen("Energias_Euler_h_0.001.txt","w");
    for(t=0;t<Nt;t++){

        tiempo=tiempo+h;
        fuerza_energia(L,particula,cutoff,fuerza,&Epot);
        Energia_cinetica(velocidad,&Ecin);
        Momento(velocidad,&p);

        E=Epot+Ecin;

        //fprintf(f3,"%lf     %lf     %lf     %lf\n",tiempo,Epot,Ecin,E);


        for(i=0;i<N;i++){
        particula[i].x=particula[i].x+velocidad[i].x*h+0.5*fuerza[i].x*h*h;
        particula[i].y=particula[i].y+velocidad[i].y*h+0.5*fuerza[i].y*h*h;
        particula[i].z=particula[i].z+velocidad[i].z*h+0.5*fuerza[i].z*h*h;

        //PERIODIC BOUNDARY CONDITIONS



        velocidad[i].x=velocidad[i].x+fuerza[i].x*h;
        velocidad[i].y=velocidad[i].y+fuerza[i].y*h;
        velocidad[i].z=velocidad[i].z+fuerza[i].z*h;

        //PERIODIC BOUNDARY CONDITIONS
        pbc(&particula[i].x,L);
        pbc(&particula[i].y,L);
        pbc(&particula[i].z,L);


        }
        fprintf(f2,"%lf    %lf     %lf\n",tiempo,E,p);

    }
    fclose(f2);
    //fclose(f3);

}



int main(){
    ini_ran(time(NULL));
    double cutoff=2.5;
    double T=100;
    double nu=0.1;
    double h=0.0001;
    double L;


    inicializar_sc(particula,&L);
    inicializar_velocidad(velocidad,T);
    Euler(particula, L,cutoff,h);
return 0;
}
