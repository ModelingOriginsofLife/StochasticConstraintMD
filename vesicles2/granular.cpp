#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

using namespace std;

#include <vector>

#include "vectormath.h"

#define DT 5e-2

#define GRAV 0
#define RADIUS 1
#define IRAD 3
#define BS 3
#define DISSIPATION 0.0
#define TYPES 4

#define HASHSIZE 10000

int RELAXATION = 1;

int canBond[TYPES*TYPES] = 
{
	0, 1, 0, 0,
	1, 1, 1, 0,
	0, 1, 1, 1,
	0, 0, 1, 0
};

int maxBonds[TYPES*TYPES] =
{
	10, 10, 10, 10,
	10, 2, 1, 10,
	10, 1, 2, 10,
	10, 10, 10, 10,	
};

float bondRatio[TYPES*TYPES] =
{
	5e-2, 5e-1, 0,    0,
	5e-1, 1e-3, 1e-5, 0,
	0,    1e-5, 1e-3, 5e-1,
	0,    0,    5e-1, 1e-2
};

float bondLengths[TYPES*TYPES] =
{
	2.5, 3,   0,     0,
	3,   3.0,   2.5,   0,
	0,   2.5,   2.2,   3,
	0,   0,     3,   2.5
};

int canRepel[TYPES*TYPES] = 
{
	0, 0, 1, 1,
	0, 0, 0, 1,
	1, 0, 0, 0,
	1, 1, 0, 0
};

int XR=120, YR=120;

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
		int NBonds[TYPES];
		
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

float getMakeReplProbability(int i, int j, float r2)
{
	//~ if (r2>pow(3,2))
	return canRepel[Parts[i].type + Parts[j].type*TYPES];
	
	return 0;	
}

float getBreakReplProbability(int i, int j)
{
	return 0.0001*canRepel[Parts[i].type + Parts[j].type*TYPES];
	
	return 0;
}

float getMakeAttrProbability(int i, int j, float r2)
{
	if (RELAXATION) return 0;
	
	if (Parts[i].NBonds[Parts[j].type]>=maxBonds[Parts[i].type+TYPES*Parts[j].type]) return 0;
	if (Parts[j].NBonds[Parts[i].type]>=maxBonds[Parts[i].type+TYPES*Parts[j].type]) return 0;
		
	if (r2<pow(bondLengths[Parts[i].type+TYPES*Parts[j].type],2))
		return canBond[Parts[i].type+Parts[j].type*TYPES];
	
	return 0;	
}

float getBreakAttrProbability(int i, int j)
{		
	return bondRatio[Parts[i].type + Parts[j].type*TYPES] * canBond[Parts[i].type+Parts[j].type*TYPES];
	
	return 0;
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
	int doCollide = 0;
	
	if (r<IRAD*IRAD) // Maximum interaction radius
	{
		//~ pair<int,int> P;
		
		//~ P.first = idx;
		//~ P.second = P2->idx;
		
		//~ bList.push_back(P);

		neighbors[P2->type]++;
		P2->neighbors[type]++;
		
		if ((!bondExists(idx, P2->idx))&&(rand()%5==0))
		{
			Bond B;
			
			B.i = idx;
			B.j = P2->idx;
			
			int pi = B.i, pj = B.j;
			
			float p = getMakeAttrProbability(pi,pj,r);
			
			//~ float p = (1.0 - pow(Parts[pi].NBonds/3.0,2))*(1.0 - pow(Parts[pj].NBonds/3.0,2)); //(Parts[pi].neighbors[1]+Parts[pj].neighbors[1]) == 1;
			
			if (rand()%10000001<10000000*p)
			{
				B.blen = bondLengths[type+TYPES*P2->type];
				B.type = 0;
				NBonds[P2->type]++; 
				P2->NBonds[type]++;
				bondGrid[B.i+B.j*Parts.size()] = 1;
				Bonds.push_back(B);
			}
			else
			{
				float p = getMakeReplProbability(pi,pj,r);
				if (rand()%10000001<10000000*p)
				{
					doCollide=1;
				}
			}
					//~ B.blen = sqrt(r);
					//~ B.type = 1;
					//~ bondGrid[B.i+B.j*Parts.size()] = 1;
					//~ Bonds.push_back(B);
				//~ }
			//~ }
		}				
	}

	float len = 2*RADIUS;
	if (doCollide)
		len = sqrt(r);
		
	if (r<(4*RADIUS*RADIUS) || doCollide)
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
		X = rMean + (len/2)*dR;
		P2->X = rMean - (len/2)*dR;
		
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

void Iterate()
{
	sortParts();
	//~ bList.clear();
	
	for (int i=0;i<Parts.size();i++)
	{
		for (int j=0;j<TYPES;j++)
			Parts[i].neighbors[j]=0;
		
		Parts[i].V.x[1] -= DT*GRAV;
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
				//~ int xm2 = (xm+XR)%XR, ym2 = (ym+YR)%YR;
				if ((xm>=0)&&(ym>=0)&&(xm<XR)&&(ym<YR))
				{
					int idx = Grid[xm+ym*XR];
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
	}
	
	
	for (int i=0;i<Bonds.size();i++)
	{
		Parts[Bonds[i].i].applyBond(&Parts[Bonds[i].j],&Bonds[i]);
	}

	// Make/break bonds randomly
	int N = 0.5*Bonds.size(); //(bList.size() + Bonds.size());
	//~ if (N<1) N=1;
	int stop=0;
	
	for (int i=0;(i<N)&&(!stop);i++)
	{
		int idx = rand()%Bonds.size();
		{
			//~ idx -= bList.size();
			
			int pi = Bonds[idx].i, pj = Bonds[idx].j;
			
			float p;
			
			if (Bonds[idx].type==0)
				p = getBreakAttrProbability(pi,pj);
			else
				p = getBreakReplProbability(pi,pj);
			
			if (rand()%10000001<10000000*p)
			{
				if (Bonds[idx].type==0)
				{
					Parts[pi].NBonds[Parts[pj].type]--; 
					Parts[pj].NBonds[Parts[pi].type]--;
				}
				
				bondGrid[Bonds[idx].i+Bonds[idx].j*Parts.size()] = 0;
				Bonds.erase(Bonds.begin()+idx);
			}
		}
		
		if (Bonds.size() == 0)
			stop = 1;
	}
}

void Init()
{
	Grid=(int*)malloc(XR*YR*sizeof(int));

	for (int i=0;i<1400*16;i++)
	{
		Part P;
		
		P.X.x[0] = XR*BS*(rand()%1000000)/1000000.0;
		P.X.x[1] = YR*BS*(rand()%1000000)/1000000.0;
		P.V.x[0] = (rand()%2000001-1000000.0)/1000000.0;
		P.V.x[1] = (rand()%2000001-1000000.0)/1000000.0;
		
		if (rand()%10 == 0)
			P.type = 1 + rand()%2;
		else P.type = 3*(rand()%4==0);
		
		P.idx = i;
		for (int j=0;j<TYPES;j++)
			P.NBonds[j] = 0;
		
		Parts.push_back(P);
	}
	
	mkBondGrid();
}

int XRes,YRes;
float *accBuf;

void addCircle(float x, float y, float rad, float r, float g, float b)
{
	int y0=y*YRes/(BS*YR), x0=x*XRes/(BS*XR), irad = rad*XRes/(BS*XR);
	if (irad<1) irad=1;
	
	for (int ym=y0-irad;ym<=y0+irad;ym++)
	{
		for (int xm=x0-irad;xm<=x0+irad;xm++)
		{
			if (pow(xm-x0,2)+pow(ym-y0,2)<=pow(irad,2))
			{
				if ((xm>=0)&&(ym>=0)&&(xm<XRes)&&(ym<YRes))
				{
					accBuf[3*(xm+ym*XRes)+0] += r;
					accBuf[3*(xm+ym*XRes)+1] += g;
					accBuf[3*(xm+ym*XRes)+2] += b;
				}
			}
		}
	}
}

void Render()
{
	for (int i=0;i<Parts.size();i++)
	{
		float r,g,b;
		
		switch (Parts[i].type)
		{
			case 0:
				r = 0.6;
				g = 0.8;
				b = 1.0;
				break;
			case 1:
				r = 0.7;
				g = 0.7;
				b = 0.3;
				break;
			case 2:
				r = 0.7;
				g = 0.3;
				b = 0.7;
				break;
			case 3:
				r = 1.0;
				g = 0.6;
				b = 0.4;
				break;
		}
		
		addCircle(Parts[i].X.x[0],Parts[i].X.x[1],1.0, r,g,b);
	}
}

int main(int argc, char **argv)
{
	int frame=0;
	Init();
	
	XRes = 200; YRes = 200; accBuf=(float*)malloc(XRes*YRes*3*sizeof(float));
	memset(accBuf,0,XRes*YRes*3*sizeof(float));
	
	RELAXATION = 1;
	
	while (1)
	{
		Iterate();
		frame++;
		if (frame%10 == 0)
			Render();
		
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
			
			f = fopen("avg.out","wb");
			for (int y=0;y<YRes;y++)
			{
				for (int x=0;x<XRes;x++)
				{
					fprintf(f,"%.6g %.6g %.6g ",accBuf[(x+y*XRes)*3],accBuf[(x+y*XRes)*3+1],accBuf[(x+y*XRes)*3+2]);
				}
				fprintf(f,"\n");
			}
			fclose(f);
			memset(accBuf,0,XRes*YRes*3*sizeof(float));
			char Str[512];
			sprintf(Str,"python -W ignore render.py frames/%.6d.png",frame/500);
			system(Str);
		}
		
		if (frame==1500) RELAXATION = 0;
	}
}
