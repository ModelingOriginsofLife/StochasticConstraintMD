#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "vectormath.h"

Vect3::Vect3()
{
    x[0]=x[1]=x[2]=0;
}

Vect3::Vect3(float xm[3])
{
    memcpy(x,xm,sizeof(float)*3);
}

Vect3::Vect3(float xm, float ym, float zm)
{
    x[0]=xm;
    x[1]=ym;
    x[2]=zm;
}

Vect3 Vect3::operator%(Vect3 B)
{
    Vect3 temp;
    int i;

    for (i=0; i<3; i++)
        temp.x[i]=x[i]*B.x[i];

    return temp;
}

Vect3 Vect3::operator+(Vect3 B)
{
    Vect3 temp;
    int i;

    for (i=0; i<3; i++)
        temp.x[i]=x[i]+B.x[i];

    return temp;
}

Vect3 Vect3::operator-(Vect3 B)
{
    Vect3 temp;
    int i;

    for (i=0; i<3; i++)
        temp.x[i]=x[i]-B.x[i];

    return temp;
}

Vect3 Vect3::operator+=(Vect3 B)
{
    (*this) = (*this) + B;
    return (*this);
}

Vect3 Vect3::operator-=(Vect3 B)
{
    (*this) = (*this) - B;
    return (*this);
}

Vect3 Vect3::operator/(Vect3 B) // Cross product!
{
    Vect3 temp;
    int i;

    for (i=0; i<3; i++)
        temp.x[i]=x[(i+1)%3]*B.x[(i+2)%3]-x[(i+2)%3]*B.x[(i+1)%3];

    return temp;
}

float Vect3::operator*(Vect3 B)
{
    float temp=0;
    int i;

    for (i=0; i<3; i++)
        temp+=x[i]*B.x[i];

    return temp;
}

Vect3 Vect3::operator/(float b)
{
    Vect3 temp;
    int i;

    for (i=0; i<3; i++)
        temp.x[i]=x[i]/b;

    return temp;
}

Vect3 operator*(float a, Vect3 B)
{
    Vect3 temp;
    int i;

    for (i=0; i<3; i++)
        temp.x[i]=B.x[i]*a;

    return temp;
}

Vect3 Vect3::operator*(float b)
{
    Vect3 temp;
    int i;

    for (i=0; i<3; i++)
        temp.x[i]=x[i]*b;

    return temp;
}

void Vect3::Normalize()
{
    float m=sqrt((*this)*(*this));

    (*this)=(*this)/m;
}

bool Vect3::Within(Vect3 A, Vect3 B)
{
    int i;
    bool out=false;

    for (i=0; (i<3)&&!out; i++)
    {
//        if (x[i]!=x[i]) return false;
//        if (x[i] == numeric_limits<float>::infinity()) return false;
        if (x[i]<=A.x[i]) out=true;
        if (x[i]>B.x[i]) out=true;
    }

    return !out;
}

void Vect3::zero()
{
    x[0]=x[1]=x[2]=0;
}

float Gaussian()
{
    float z1,z2;

    z1=(rand()%10000000+1)/10000000.0;
    z2=(rand()%10000000+1)/10000000.0;

    return sqrt(-2*log(z1))*cos(2*M_PI*z2);
}

void Vect3::GaussRandom(float stdev)
{
    for (int i=0; i<3; i++) x[i]=Gaussian()*stdev;
}

void Vect3::BoxRandom(Vect3 A, Vect3 B)
{
    for (int i=0; i<3; i++) x[i]=A.x[i]+(B.x[i]-A.x[i])*(rand()%10000000)/10000000.0;
}

// 2D routines

Vect2::Vect2()
{
    x[0]=x[1]=0;
}

Vect2::Vect2(float xm[2])
{
    memcpy(x,xm,sizeof(float)*2);
}

Vect2::Vect2(float xm, float ym)
{
    x[0]=xm;
    x[1]=ym;
}

Vect2 Vect2::operator%(Vect2 B)
{
    Vect2 temp;
    int i;

    for (i=0; i<2; i++)
        temp.x[i]=x[i]*B.x[i];

    return temp;
}

Vect2 Vect2::operator+(Vect2 B)
{
    Vect2 temp;
    int i;

    for (i=0; i<2; i++)
        temp.x[i]=x[i]+B.x[i];

    return temp;
}

Vect2 Vect2::operator-(Vect2 B)
{
    Vect2 temp;
    int i;

    for (i=0; i<2; i++)
        temp.x[i]=x[i]-B.x[i];

    return temp;
}

Vect2 Vect2::operator+=(Vect2 B)
{
    (*this) = (*this) + B;
    return (*this);
}

Vect2 Vect2::operator-=(Vect2 B)
{
    (*this) = (*this) - B;
    return (*this);
}

float Vect2::operator/(Vect2 B) // Cross product!
{
    int i;

    return x[0]*B.x[1]-x[1]*B.x[0];
}

float Vect2::operator*(Vect2 B)
{
    float temp=0;
    int i;

    for (i=0; i<2; i++)
        temp+=x[i]*B.x[i];

    return temp;
}

Vect2 Vect2::operator/(float b)
{
    Vect2 temp;
    int i;

    for (i=0; i<2; i++)
        temp.x[i]=x[i]/b;

    return temp;
}

Vect2 operator*(float a, Vect2 B)
{
    Vect2 temp;
    int i;

    for (i=0; i<2; i++)
        temp.x[i]=B.x[i]*a;

    return temp;
}

Vect2 Vect2::operator*(float b)
{
    Vect2 temp;
    int i;

    for (i=0; i<2; i++)
        temp.x[i]=x[i]*b;

    return temp;
}

void Vect2::Normalize()
{
    float m=sqrt((*this)*(*this));

    (*this)=(*this)/m;
}

bool Vect2::Within(Vect2 A, Vect2 B)
{
    int i;
    bool out=false;

    for (i=0; (i<2)&&!out; i++)
    {
//        if (x[i]!=x[i]) return false;
//        if (x[i] == numeric_limits<float>::infinity()) return false;
        if (x[i]<=A.x[i]) out=true;
        if (x[i]>B.x[i]) out=true;
    }

    return !out;
}

void Vect2::zero()
{
    x[0]=x[1]=x[2]=0;
}


void Vect2::GaussRandom(float stdev)
{
    for (int i=0; i<2; i++) x[i]=Gaussian()*stdev;
}

void Vect2::BoxRandom(Vect2 A, Vect2 B)
{
    for (int i=0; i<2; i++) x[i]=A.x[i]+(B.x[i]-A.x[i])*(rand()%10000000)/10000000.0;
}

Mat33::Mat33()
{
    int i,j;

    for (j=0; j<3; j++)
        for (i=0; i<3; i++)
            x[i][j]=0;
}

Mat33::Mat33(Vect3 W)
{
    x[0][0]=0;       x[1][0]=W.x[2];  x[2][0]=-W.x[1];
    x[0][1]=-W.x[2]; x[1][1]=0;       x[2][1]=W.x[0];
    x[0][2]=W.x[1];  x[1][2]=-W.x[0]; x[2][2]=0;
}

Mat33 Mat33::operator+(Mat33 B)
{
    Mat33 temp;
    int i,j;

    for (j=0; j<3; j++)
        for (i=0; i<3; i++)
            temp.x[i][j]=x[i][j]+B.x[i][j];

    return temp;
}

Mat33 Mat33::operator-(Mat33 B)
{
    Mat33 temp;
    int i,j;

    for (j=0; j<3; j++)
        for (i=0; i<3; i++)
            temp.x[i][j]=x[i][j]-B.x[i][j];

    return temp;
}

Mat33 Mat33::operator+=(Mat33 B)
{
    (*this) = (*this) + B;
    return (*this);
}

Mat33 Mat33::operator-=(Mat33 B)
{
    (*this) = (*this) - B;
    return (*this);
}

Vect3 Mat33::operator*(Vect3 B)
{
    Vect3 temp(0,0,0);
    int i,j;

    for (j=0; j<3; j++)
        for (i=0; i<3; i++)
            temp.x[j]+=x[i][j]*B.x[i];

    return temp;
}

Mat33 Mat33::operator*(Mat33 B)
{
    Mat33 temp;
    int i,j,k;

    for (j=0; j<3; j++)
        for (i=0; i<3; i++)
            for (k=0; k<3; k++)
                temp.x[i][j]+=x[k][j]*B.x[i][k];

    return temp;
}

Mat33 Mat33::operator/(float b)
{
    Mat33 temp;
    int i,j;

    for (j=0; j<3; j++)
        for (i=0; i<3; i++)
            temp.x[i][j]=x[i][j]/b;

    return temp;
}

Mat33 operator*(float a, Mat33 B)
{
    Mat33 temp;
    int i,j;

    for (j=0; j<3; j++)
        for (i=0; i<3; i++)
            temp.x[i][j]=B.x[i][j]*a;

    return temp;
}

Mat33 Mat33::operator*(float b)
{
    Mat33 temp;
    int i,j;

    for (j=0; j<3; j++)
        for (i=0; i<3; i++)
            temp.x[i][j]=x[i][j]*b;

    return temp;
}

void Mat33::zero()
{
    for (int j=0; j<3; j++)
        for (int i=0; i<3; i++)
            x[i][j]=0;
}

void Mat33::ident()
{
    for (int j=0; j<3; j++)
        for (int i=0; i<3; i++)
            x[i][j]=(i==j);
}
