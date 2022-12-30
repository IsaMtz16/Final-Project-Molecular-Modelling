#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*N:number of particles
 Nt:number of time steps
*/

#define N 125
#define Nt 100000

extern float Random(void);
#define NormRANu (2.3283063671E-10)
extern void ini_ran(int SEMILLA);
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran, ig1, ig2, ig3;

#define PI 3.14159265

/*First we create a new type of variable*/
typedef struct{
    double x, y, z;
} Vector;

/*We declare globally two arrays of type Vector of
dimension N for the position and the velocity*/
Vector particula[N]={0};
Vector velocidad[N]={0};

//Random number generator
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
/*Function for initialize the position of the
particles in simple cubic lattice*/

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
/*Function to initialise the velocity with a bimodal distribution */
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
/*Function to generate random numbers with a gaussian distribution*/
void Box_Muller (double *g1, double *g2,double sigma)
{
    double d1, d2;
    d1=sigma*sqrt(-2.0*log(Random()));
    d2=2.0*PI*Random();
    *g1=-d1*cos(d2);
    *g2=-d1*sin(d2);
}

/*Function of the Andersen thermostat*/
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

/*Function to implement periodic boundary conditions
We have a box of size L*L centered on zero
After updating the position of the particles we check if any of the components
of the position vector
*/

void pbc(double *x,double L){
    if(*x>(L/2)){
        *x=*x-L;
    }
    if(*x<-(L/2)){
        *x=*x+L;
    }

}
/*Function for initialising all the components of the arrays of type Vector of dimension N
to zero*/
void vector_0(Vector *fuerza){
    int i;
    for(i=0;i<N;i++){
        fuerza[i].x=0;
        fuerza[i].y=0;
        fuerza[i].z=0;

    }
}
/*Function for computing the force that suffers each particle and the total energy considering
the system in a Lennard-Jones potential*/
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
/*Function for obtaining the kinetic energy of the system at a certain
time*/
void Energia_cinetica(Vector velocidad[],double *Ecin){
    *Ecin=0;
    int i;
    for(i=0;i<N;i++){
        *Ecin=*Ecin+velocidad[i].x*velocidad[i].x+velocidad[i].y*velocidad[i].y+velocidad[i].z*velocidad[i].z;
    }
    *Ecin=*Ecin/2;
}
/*Function for calculating the total moment of the system*/
void Momento(Vector velocidad[],double *p){
    int i;
    *p=0;
    for(i=0;i<N;i++){
        *p=*p+sqrt(velocidad[i].x*velocidad[i].x+velocidad[i].y*velocidad[i].y+velocidad[i].z*velocidad[i].z);
    }
}
/*Function for computing the temperature of the system*/
void Temperatura(double Ecin,double *T_inst){
        *T_inst=2*Ecin/(3*(double)N-3);
}
/*Function to calculate the module of a vector*/
void modulo(Vector velocidad[],int i,double *v){
    *v=sqrt(velocidad[i].x*velocidad[i].x+velocidad[i].y*velocidad[i].y+velocidad[i].z*velocidad[i].z);
}

/*Implementation of the algorithm of velocity Verlet*/
void velocity_Verlet(double h,double L,Vector *particula,double cutoff,double T,double nu){
    int i,j,k,m;
    double t=0;
    double T_inst;
    double E, Epot, Ecin;
    double p;//momento
    double media_T=0;
    E=0;
    double v;//modulo de la velocidad


    //ABRIMOS UN FICHERO PARA ESCRIBIR LAS TRAYECTORIAS
    /*FILE*f;
    f=fopen("Trayectoria_Verlet_h_0.1.txt","w");*/

    //ABRIMOS UN FICHERO PARA ESCRIBIR LAS ENERGÍAS
    //FILE*f1;
    //f1=fopen("Energias_Verlet_h_0.001.txt","w");
    Vector fuerza[N]={0};

    /*FILE*f2;
    f2=fopen("Energia_momento_t_Verlet_h_0.01.txt","w");*/

    /*FILE*f3;
    f3=fopen("Temperatura_t.txt","w");*/

    /*FILE*f4;
    f4=fopen("Distr_v_inicial2.txt","w");*/
    for(i=0;i<Nt;i++){

        t=t+h;
        fuerza_energia(L,particula,cutoff,fuerza,&Epot);
        Energia_cinetica(velocidad,&Ecin);
        E=Epot+Ecin;
        Momento(velocidad,&p);
        /*if((i%1)==0){
            fprintf(f1,"%lf     %lf     %lf     %lf\n",t,Epot,Ecin,E);
        }*/


        //fprintf(f,"%d\n\n",N);
        //vector_0(velocidad);

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




            //fprintf(f,"S    %lf     %lf     %lf\n",particula[j].x,particula[j].y,particula[j].z);
        }
        fuerza_energia(L,particula,cutoff,fuerza,&Epot);
        for(k=0;k<N;k++){

            velocidad[k].x=velocidad[k].x+0.5*fuerza[k].x*h;
            velocidad[k].y=velocidad[k].y+0.5*fuerza[k].y*h;
            velocidad[k].z=velocidad[k].z+0.5*fuerza[k].z*h;

            if(i==0){
                modulo(velocidad,k,&v);
                //fprintf(f4,"%lf\n",v);
            }

            //Anderson_thermostat(velocidad,T,nu);
        }
        Energia_cinetica(velocidad,&Ecin);
        E=Epot+Ecin;
        Temperatura(Ecin,&T_inst);
        media_T=media_T+T_inst;
        Momento(velocidad,&p);
        //printf("%d\n",i);

        if(i%10==0){
            //fprintf(f2,"%lf     %lf     %lf\n",t,E,p);
            //fprintf(f3,"%lf     %lf\n",t,T_inst);
        }







    }
    media_T=media_T/(double)Nt;
    FILE*f6;
    f6=fopen("media_T.txt","w");
    fprintf(f6,"%lf \n",media_T);
    fclose(f6);
    /*FILE*f5;
    f5=fopen("Distr_v_final2.txt","w");
    for(m=0;m<N;m++){
        modulo(velocidad,m,&v);
        fprintf(f5,"%lf\n",v);
    }*/


    //fclose(f2);
    //fclose(f3);
    //fclose(f4);
    //fclose(f5);

    //fclose(f);
    //fclose(f1);
}

int main(){
    ini_ran(time(NULL));
    double cutoff=2.5;
    double T=100;
    double nu=0.1;
    double h=0.001;
    double L;


    inicializar_sc(particula,&L);
    inicializar_velocidad(velocidad,T);
    velocity_Verlet(h,L,particula,cutoff,T,nu);
return 0;
}
