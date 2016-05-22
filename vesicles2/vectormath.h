class Vect3
{
public:
    float x[3];

    Vect3();
    Vect3(float xm[3]);
    Vect3(float xm, float ym, float zm);
    Vect3 operator +(Vect3 B);
    Vect3 operator -(Vect3 B);
    Vect3 operator+=(Vect3 B);
    Vect3 operator-=(Vect3 B);
    Vect3 operator *(float b); // scale
    Vect3 operator /(float b); // scale
    Vect3 operator /(Vect3 B); // cross product
    float operator *(Vect3 B); // dot product
    Vect3 operator %(Vect3 B); // piecewise product
    void Normalize();
    bool Within(Vect3 A, Vect3 B); // Bounding box check
    void zero();
    void GaussRandom(float stdev);
    void BoxRandom(Vect3 A, Vect3 B);
};

extern Vect3 operator*(float a, Vect3 B);
extern float Gaussian();

class Vect2
{
public:
    float x[2];

    Vect2();
    Vect2(float xm[2]);
    Vect2(float xm, float ym);
    Vect2 operator +(Vect2 B);
    Vect2 operator -(Vect2 B);
    Vect2 operator+=(Vect2 B);
    Vect2 operator-=(Vect2 B);
    Vect2 operator *(float b); // scale
    Vect2 operator /(float b); // scale
    float operator /(Vect2 B); // cross product
    float operator *(Vect2 B); // dot product
    Vect2 operator %(Vect2 B); // piecewise product
    void Normalize();
    bool Within(Vect2 A, Vect2 B); // Bounding box check
    void zero();
    void GaussRandom(float stdev);
    void BoxRandom(Vect2 A, Vect2 B);
};

extern Vect2 operator*(float a, Vect2 B);

class Mat33
{
public:
    float x[3][3]; // components

    Mat33();
    Mat33(Vect3 W); // Rotation matrix

    Mat33 operator +(Mat33 B);
    Mat33 operator -(Mat33 B);
    Mat33 operator+=(Mat33 B);
    Mat33 operator-=(Mat33 B);
    Mat33 operator *(float b); // scale
    Mat33 operator /(float b); // scale
    Mat33 operator *(Mat33 B); // matrix product
    Vect3 operator *(Vect3 B); // matrix times vector

    void zero();
    void ident();
};

extern Mat33 operator*(float a, Mat33 B);
