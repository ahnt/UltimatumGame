//
//  main.cpp
//  CrowdSourcing
//
//  Created by Arend Hintze on 2/6/13.
//  Copyright (c) 2013 MPI-BERLIN. All rights reserved.
//

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <time.h>

#define randDouble ((double)rand()/(double)RAND_MAX)

int globalUpdate=0;
using namespace std;
double u=0.001; //the mutation rate
double w=10.0; //selection strength
int popSize=100;
double limit;

class tAgent{
public:
    tAgent *ancestor;
    int nrPointingAtMe;
    int born;
    double genome[3][2];
    tAgent();
    ~tAgent();
    void setupRand(int startCondition);
    void inherit(tAgent *from);
    void LOD(FILE *F);
};

vector<tAgent*> population;
vector<vector<double> > payoffs;
vector<double> sumOfPayoffs;
vector<double> fitness;

void recalculateEverything(void);
void showPayoffs(void);
void popCheck(void);
void recalculateSingle(int who);

int main(int argc, const char * argv[])
{
    int i,j,g,k,x,y,nx,ny,deadGuy,newGuy;
    double maxFit;
    srand((int)time(NULL));
//    srand((int)getpid());
	u=atof(argv[2]);
	w=atof(argv[3]);
    limit=atof(argv[4]);
	popSize=atoi(argv[5]);
	printf("U:%f W:%f popSize\n",u,w,popSize);
    population.clear();
	payoffs.resize(popSize);
	sumOfPayoffs.resize(popSize);
	fitness.resize(popSize);
    for(i=0;i<popSize;i++){
		payoffs[i].resize(popSize);
        tAgent *A=new tAgent;
        A->setupRand(-1);
        population.push_back(A);
    }
    recalculateEverything();
    for(globalUpdate=1;globalUpdate<5000000;globalUpdate++){
        //showPayoffs();
        maxFit=0.0;
        for(i=0;i<popSize;i++){
            fitness[i]=exp(w*(sumOfPayoffs[i]/(double)(popSize-1)));
            if(maxFit<fitness[i])
                maxFit=fitness[i];
        }
        do{
            newGuy=rand()%popSize;
        }while(randDouble>(fitness[newGuy]/maxFit));
        do{
            deadGuy=rand()%popSize;
        }while(deadGuy==newGuy);
        population[deadGuy]->nrPointingAtMe--;
        if(population[deadGuy]->nrPointingAtMe==0)
            delete population[deadGuy];
        population[deadGuy]=new tAgent;
        population[deadGuy]->inherit(population[newGuy]);
        recalculateSingle(deadGuy);
        if((globalUpdate&16383)==0){
            cout<<globalUpdate<<" ";
            popCheck();
            cout<<endl;
        }
    }
    FILE *F=fopen(argv[1],"w+t");
    population[0]->LOD(F);
    fclose(F);
    return 0;
}

tAgent::tAgent(){
    ancestor=NULL;
    nrPointingAtMe=1;
    born=globalUpdate;
}

tAgent::~tAgent(){
    if(ancestor!=NULL){
        ancestor->nrPointingAtMe--;
        if(ancestor->nrPointingAtMe==0)
            delete ancestor;
    }
}

void tAgent::setupRand(int startCondition){
    int i,j;
	if(startCondition==-1){
		for(i=0;i<3;i++)
			for(j=0;j<2;j++)
				genome[i][j]=randDouble;
		genome[2][1]=limit+(randDouble*(1.0-limit));
	} else {
		int z=0;
		for(i=0;i<3;i++)
			for(j=0;j<2;j++){
				genome[i][j]=(double)((startCondition>>z)&1);
				z++;
			}
	}
}

void tAgent::inherit(tAgent *from){
    int i,j;
    from->nrPointingAtMe++;
    ancestor=from;
    if(randDouble<u){
        for(i=0;i<3;i++)
            for(j=0;j<2;j++){
                genome[i][j]=randDouble;
			}
        genome[2][1]=limit+(randDouble*(1.0-limit));
    }
    else{
        for(i=0;i<3;i++)
            for(j=0;j<2;j++)
                genome[i][j]=from->genome[i][j];
    }
}

void tAgent::LOD(FILE *F){
    if(ancestor!=NULL)
        ancestor->LOD(F);
    else
        fprintf(F,"born,ps,qs,po,qo,id,th\n");
    fprintf(F,"%i,%f,%f,%f,%f,%f,%f\n",born
            ,genome[0][0],genome[0][1]
            ,genome[1][0],genome[1][1]
            ,genome[2][0],genome[2][1]);
}

// *** regular functions

void recalculateEverything(void){
    int i,j;
    for(i=0;i<popSize;i++){
        fitness[i]=0.0;
        sumOfPayoffs[i]=0.0;
        for(j=0;j<popSize;j++)
            payoffs[i][j]=0.0;
    }
    for(i=0;i<popSize;i++)
        for(j=i+1;j<popSize;j++){
            int s1,s2;
            if(fabs(population[i]->genome[2][0]-population[j]->genome[2][0])<population[i]->genome[2][1])
                s1=0;
            else
                s1=1;
            if(fabs(population[i]->genome[2][0]-population[j]->genome[2][0])<population[j]->genome[2][1])
                s2=0;
            else
                s2=1;
            if(population[i]->genome[s1][0]>population[j]->genome[s2][1]){
                payoffs[i][j]+=1.0-population[i]->genome[s1][0];
                payoffs[j][i]+=population[i]->genome[s1][0];
            }
            //j->i
            if(population[j]->genome[s2][0]>population[i]->genome[s1][1]){
                payoffs[j][i]+=1.0-population[j]->genome[s2][0];
                payoffs[i][j]+=population[j]->genome[s2][0];
            }
        }
    for(i=0;i<popSize;i++){
        for(j=0;j<popSize;j++)
            if(i!=j)
                sumOfPayoffs[i]+=payoffs[i][j];
    }
}

void recalculateSingle(int who){
    int i,j;
    sumOfPayoffs[who]=0.0;
    for(i=0;i<popSize;i++){
        sumOfPayoffs[i]-=payoffs[i][who];
        payoffs[i][who]=0.0;
        payoffs[who][i]=0.0;
    }
    for(i=0;i<popSize;i++)
        if(i!=who){
            int s1,s2;
            if(fabs(population[i]->genome[2][0]-population[who]->genome[2][0])<population[i]->genome[2][1])
                s1=0;
            else
                s1=1;
            if(fabs(population[i]->genome[2][0]-population[who]->genome[2][0])<population[who]->genome[2][1])
                s2=0;
            else
                s2=1;
            if(population[i]->genome[s1][0]>population[who]->genome[s2][1]){
                payoffs[i][who]+=1.0-population[i]->genome[s1][0];
                payoffs[who][i]+=population[i]->genome[s1][0];
            }
            //j->i
            if(population[who]->genome[s2][0]>population[i]->genome[s1][1]){
                payoffs[who][i]+=1.0-population[who]->genome[s2][0];
                payoffs[i][who]+=population[who]->genome[s2][0];
            }        
        }
    for(i=0;i<popSize;i++)
        if(i!=who){
            sumOfPayoffs[i]+=payoffs[i][who];
            sumOfPayoffs[who]+=payoffs[who][i];
        }
}
void showPayoffs(void){
    int i,j;
    for(j=0;j<popSize;j++){
        for(i=0;i<popSize;i++)
            printf("%0.02f ",payoffs[i][j]);
        printf("\n");
    }
}

void popCheck(void){
    int i,j=0;
    double meanP=0.0,meanF=0.0,maxFit=0.0;
    for(i=0;i<popSize;i++){
        meanP+=sumOfPayoffs[i]/(double)(popSize-1);
        meanF+=fitness[i];
        if(fitness[i]>maxFit){
            j=i;
            maxFit=fitness[i];
        }
    }
    cout<<(meanP/(double)(popSize-1))<<" "<<(meanF/(double)(popSize-1))<<" "<<population[j]->genome[0][0]<<" "<<population[j]->genome[0][1];
}

