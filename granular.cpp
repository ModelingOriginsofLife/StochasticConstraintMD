#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

using namespace std;

#include <vector>

#include "vectormath.h"

#define DT 1e-2

#define RADIUS 1
#define BOND 3
#define BS 3
#define DISSIPATION 0.0
#define TYPES 4

vector< pair<int,int> > bList;

int XR=40, YR=40;

class Part;

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
		void applyBond(Part *P2);
		void Bound();
};

class Bond
{
	public:
		int i,j;
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
	return bondGrid[i+j*Parts.size()];
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
	
	if (r<BOND*BOND)
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

void Part::applyBond(Part *P2)
{
	Vect2 dR = X-P2->X;
	float r = dR*dR;
	
	if (r>BOND*BOND)
	{
		if (r<1e-3) r=1e-3;
		r=sqrt(r);

		dR = dR/r;
		
		Vect2 VCM=0.5*(V+P2->V);
		Vect2 vIn = V-VCM;
		
		float F=-(2-DISSIPATION)*vIn*dR;
		if (F>0) F=0;
		
		Vect2 vOut = vIn+F*dR;
		
		V = VCM + vOut;
		P2->V = VCM - vOut;

		Vect2 rMean = 0.5*(X+P2->X);
		X = rMean + (BOND/2.0)*dR;
		P2->X = rMean - (BOND/2.0)*dR;
		
		Bound();
		P2->Bound();
	}
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
		Parts[Bonds[i].i].applyBond(&Parts[Bonds[i].j]);
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
				//~ Vect2 dR = Parts[B.i].X - Parts[B.j].X;
				//~ float r = dR*dR;				
				//~ float p = 0.008*r*r;
				
				float p = (1.0 - pow(Parts[pi].NBonds/3.0,2))*(1.0 - pow(Parts[pj].NBonds/3.0,2)); //(Parts[pi].neighbors[1]+Parts[pj].neighbors[1]) == 1;
				
				if (rand()%1000001<1000000*p)
				{
					Parts[pi].NBonds++; Parts[pj].NBonds++;
					bondGrid[B.i+B.j*Parts.size()] = 1;
					Bonds.push_back(B);
				}
			}
			
			bList.erase(bList.begin()+idx);
		}
		else
		{
			idx -= bList.size();
			
			int pi = Bonds[idx].i, pj = Bonds[idx].j;
			//~ Vect2 dR = Parts[Bonds[idx].i].X - Parts[Bonds[idx].j].X;
			//~ float r = dR*dR;			
			//~ float p = 0.0008*r*r;
			
			float p = 0.1 * (1.0 - pow( (Parts[pi].NBonds-1)/3.0,2))*(1.0 - pow( (Parts[pj].NBonds-1)/3.0,2)); //= 0.05*( (Parts[pi].neighbors[1]+Parts[pj].neighbors[1]) == 1);
			
			if (rand()%1000001<1000000*p)
			{
				Parts[pi].NBonds--; Parts[pj].NBonds--;
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

	for (int i=0;i<1000;i++)
	{
		Part P;
		
		P.X.x[0] = XR*BS*(rand()%1000000)/1000000.0;
		P.X.x[1] = YR*BS*(rand()%1000000)/1000000.0;
		P.V.x[0] = (rand()%2000001-1000000.0)/1000000.0;
		P.V.x[1] = (rand()%2000001-1000000.0)/1000000.0;
		P.type = rand()%TYPES;
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
		
		if (frame%1000 == 0)
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
				int j = Bonds[i].i;
				int k = Bonds[i].j;
				fprintf(f,"%.6g %.6g %.6g %.6g\n",Parts[j].X.x[0],Parts[j].X.x[1],Parts[k].X.x[0],Parts[k].X.x[1]);
			}
			fclose(f);
			
			char Str[512];
			sprintf(Str,"python -W ignore render.py frames/%.6d.png &> /dev/null",frame/1000);
			system(Str);
		}
	}
}
