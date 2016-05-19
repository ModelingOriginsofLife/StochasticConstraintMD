#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

using namespace std;

#include <vector>

#include "vectormath.h"

#define DT 1e-2

#define RADIUS 1
#define IRAD 5.0
#define BS 5.0
#define DISSIPATION 0.0
#define TYPES 3

float Interaction[TYPES*TYPES] = {
	0.0001, 0.00001, 0.00005,
	0.00001, 1.0,     0.0001,
	0.00005, 0.0001,  0.1
};

int bondPref[TYPES] = { 1, 2, 3 };

vector< pair<int,int> > bList;

int XR=20, YR=20;

class Part;

class Bond
{
	public:
		int i,j;
		int type;
		float blen;
};

class Part
{
	public:
		Vect2 X,V;		
		int type;
		int neighbors[TYPES];
		int next;
		int idx;
		int NBonds;
		
		void getGridPos(int &x, int &y);
		void tryCollide(Part *P2);
		void applyBond(Part *P2, Bond *B);
		void Bound();
};

vector<Part> Parts;
vector<Bond> Bonds;
int *Grid;
char *bondGrid;

void mkBondGrid()
{
	bondGrid=(char*)malloc(Parts.size()*Parts.size());
	memset(bondGrid,0,Parts.size()*Parts.size());
}

int bondExists(int i, int j)
{
	return (bondGrid[i+j*Parts.size()]!=0);
}

void Part::getGridPos(int &x, int &y)
{
	x = X.x[0]/BS;
	y = X.x[1]/BS;
	
	if (x<0)
		printf("%d %.6g\n",x,X.x[0]);

	if (y<0)
		printf("%d %.6g\n",y,X.x[1]);
}

void sortParts()
{
	memset(Grid,0,XR*YR*sizeof(int));
	
	for (int i=0;i<Parts.size();i++)
	{
		int x,y;
		
		Parts[i].getGridPos(x,y);
		
		if (!Grid[x+y*XR])
		{
			Parts[i].next = 0;
		}
		else
		{
			Parts[i].next = Grid[x+y*XR];
		}
		Grid[x+y*XR] = i+1;
	}
}

void Part::Bound()
{
	for (int j=0;j<2;j++)
	{
		if (X.x[j]<RADIUS)
		{
			X.x[j]=RADIUS;
			if (V.x[j]<0)
				V.x[j]*=-1;
		}

		if (X.x[j]>BS*XR-RADIUS)
		{
			X.x[j]=BS*XR-RADIUS;
			if (V.x[j]>0)
				V.x[j]*=-1;
		}
	}
}

void Part::tryCollide(Part *P2)
{
	Vect2 dR = X-P2->X;
	float r = dR*dR;
	
	if (r<IRAD*IRAD) // Maximum interaction radius
	{
		pair<int,int> P;
		
		P.first = idx;
		P.second = P2->idx;
		
		bList.push_back(P);

		neighbors[P2->type]++;
		P2->neighbors[type]++;
	}

	if (r<(4*RADIUS*RADIUS))
	{
		if (r<1e-3) r=1e-3;
		r=sqrt(r);
		
		dR = dR/r;
				
		Vect2 VCM=0.5*(V+P2->V);
		Vect2 vIn = V-VCM;
		
		float F=-(2-DISSIPATION)*vIn*dR;
		if (F<0) F=0;
		
		Vect2 vOut = vIn+F*dR;
		
		V = VCM + vOut;
		P2->V = VCM - vOut;

		Vect2 rMean = 0.5*(X+P2->X);
		X = rMean + RADIUS*dR;
		P2->X = rMean - RADIUS*dR;
		
		Bound();
		P2->Bound();
	}	
}

void Part::applyBond(Part *P2, Bond *B)
{
	Vect2 dR = X-P2->X;
	float r = dR*dR;
	float del = r-pow(B->blen,2);
	
	if (((B->type==0)&&(del>0) )||( (B->type==1)&&(del<0)))
	{
		if (r<1e-3) r=1e-3;
		r=sqrt(r);

		dR = dR/r;
		
		Vect2 VCM=0.5*(V+P2->V);
		Vect2 vIn = V-VCM;
		
		float F=-(2-DISSIPATION)*vIn*dR;
		
		if ((B->type == 0)&&(F>0))
			F=0;
		if ((B->type == 1)&&(F<0))
			F=0;
		
		Vect2 vOut = vIn+F*dR;
		
		V = VCM + vOut;
		P2->V = VCM - vOut;

		Vect2 rMean = 0.5*(X+P2->X);
		X = rMean + (B->blen/2.0)*dR;
		P2->X = rMean - (B->blen/2.0)*dR;
		
		Bound();
		P2->Bound();
	}
}

float getMakeReplProbability(int i, int j)
{
	Vect2 dR = Parts[i].X-Parts[j].X;
	float r = dR*dR;
	
	if (Parts[i].NBonds != bondPref[Parts[i].type]) return 0;
	if (Parts[j].NBonds != bondPref[Parts[j].type]) return 0;
	
	if (r>pow(4.2,2))
		return 1;
	
	return 0;	
}

float getBreakReplProbability(int i, int j)
{
	Vect2 dR = Parts[i].X-Parts[j].X;
	float r = dR*dR;
	
	if (r>pow(5,2)) return 1;
	else if (r>pow(4.2,2)) return 0.005;
	//~ else 
		//~ return 0.01*getMakeReplProbability(i,j);
}

float getMakeAttrProbability(int i, int j)
{
	Vect2 dR = Parts[i].X-Parts[j].X;
	float r = dR*dR;
	
	if (Parts[i].NBonds+1 > bondPref[Parts[i].type]) return 0;
	if (Parts[j].NBonds+1 > bondPref[Parts[j].type]) return 0;
	
	if (r<pow(2.4,2))
		return 1.0;
	
	return 0;	
}

float getBreakAttrProbability(int i, int j)
{
	Vect2 dR = Parts[i].X-Parts[j].X;
	float r = dR*dR;

	if (Parts[i].NBonds-1 < bondPref[Parts[i].type]-1) return 0;
	if (Parts[j].NBonds-1 < bondPref[Parts[j].type]-1) return 0;

	float p0 = Interaction[Parts[i].type+TYPES*Parts[i].type];

	if (r<pow(2.4,2))
		return p0 * 0.05;
	
	return 0;
}

void Iterate()
{
	sortParts();
	bList.clear();
	
	for (int i=0;i<Parts.size();i++)
	{
		for (int j=0;j<TYPES;j++)
			Parts[i].neighbors[j]=0;
		
		Parts[i].X = Parts[i].X + DT*Parts[i].V;
		Parts[i].Bound();
	}

	for (int i=0;i<Parts.size();i++)
	{
		int x,y,xm,ym;
		
		Parts[i].getGridPos(x,y);
		
		for (int ym=y-1;ym<=y+1;ym++)
			for (int xm=x-1;xm<=x+1;xm++)
			{
				int xm2 = (xm+XR)%XR, ym2 = (ym+YR)%YR;
				
				int idx = Grid[xm2+ym2*XR];
				while (idx>0)
				{
					if (idx-1<i)
					{
						Parts[i].tryCollide(&Parts[idx-1]);
					}
					idx = Parts[idx-1].next;
				}
			}
	}
	
	
	for (int i=0;i<Bonds.size();i++)
	{
		Parts[Bonds[i].i].applyBond(&Parts[Bonds[i].j],&Bonds[i]);
	}

	// Make/break bonds randomly
	int N = bList.size() + Bonds.size();
	int stop=0;
	
	for (int i=0;(i<N)&&(!stop);i++)
	{
		int idx = rand()%(bList.size()+Bonds.size());
		
		if (idx<bList.size()) // Try making a bond
		{
			if (!bondExists(bList[idx].first,bList[idx].second))
			{
				Bond B;
				
				B.i = bList[idx].first;
				B.j = bList[idx].second;
				
				int pi = B.i, pj = B.j;
				
				float p = getMakeAttrProbability(pi,pj);
				
				//~ float p = (1.0 - pow(Parts[pi].NBonds/3.0,2))*(1.0 - pow(Parts[pj].NBonds/3.0,2)); //(Parts[pi].neighbors[1]+Parts[pj].neighbors[1]) == 1;
				
				if (rand()%1000001<1000000*p)
				{
					B.blen = 2.4;
					B.type = 0;
					Parts[pi].NBonds++; Parts[pj].NBonds++;
					bondGrid[B.i+B.j*Parts.size()] = 1;
					Bonds.push_back(B);
				}
				else
				{
					float p = getMakeReplProbability(pi,pj);
					if (rand()%1000001<1000000*p)
					{
						B.blen = 4.2;
						B.type = 1;
						//~ Parts[pi].NBonds++; Parts[pj].NBonds++;
						bondGrid[B.i+B.j*Parts.size()] = 1;
						Bonds.push_back(B);
					}
				}
			}
			
			bList.erase(bList.begin()+idx);
		}
		else
		{
			idx -= bList.size();
			
			int pi = Bonds[idx].i, pj = Bonds[idx].j;
			
			float p;
			
			if (Bonds[idx].type==0)
				p = getBreakAttrProbability(pi,pj);
			else
				p = getBreakReplProbability(pi,pj);
			
			if (rand()%1000001<1000000*p)
			{
				if (Bonds[idx].type==0)
				{
					Parts[pi].NBonds--; Parts[pj].NBonds--;
				}
				
				bondGrid[Bonds[idx].i+Bonds[idx].j*Parts.size()] = 0;
				Bonds.erase(Bonds.begin()+idx);
			}
		}
		
		if (bList.size()+Bonds.size() == 0)
			stop = 1;
	}
}

void Init()
{
	Grid=(int*)malloc(XR*YR*sizeof(int));

	for (int i=0;i<400;i++)
	{
		Part P;
		
		P.X.x[0] = XR*BS*(rand()%1000000)/1000000.0;
		P.X.x[1] = YR*BS*(rand()%1000000)/1000000.0;
		P.V.x[0] = (rand()%2000001-1000000.0)/1000000.0;
		P.V.x[1] = (rand()%2000001-1000000.0)/1000000.0;
		P.type = rand()%TYPES;
		if (rand()%4 != 0)
			P.type = 0;
			
		P.idx = i;
		P.NBonds = 0;
		
		Parts.push_back(P);
	}
	
	mkBondGrid();
}

int main(int argc, char **argv)
{
	int frame=0;
	Init();
	
	while (1)
	{
		Iterate();
		frame++;
		
		if (frame%500 == 0)
		{
			FILE *f=fopen("parts.out","wb");
			for (int i=0;i<Parts.size();i++)
			{
				fprintf(f,"%.6g %.6g %d\n",Parts[i].X.x[0],Parts[i].X.x[1],Parts[i].type);
			}
			fclose(f);
			
			f = fopen("bonds.out","wb");
			for (int i=0;i<Bonds.size();i++)
			{
				if (Bonds[i].type==0)
				{
					int j = Bonds[i].i;
					int k = Bonds[i].j;
					fprintf(f,"%.6g %.6g %.6g %.6g\n",Parts[j].X.x[0],Parts[j].X.x[1],Parts[k].X.x[0],Parts[k].X.x[1]);
				}
			}
			fclose(f);
			
			char Str[512];
			sprintf(Str,"python -W ignore render.py frames/%.6d.png",frame/500);
			system(Str);
		}
	}
}
