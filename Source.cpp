
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////  10/01/2019                    Nagaraj Raparthi - VIZA 689 - Homework 2                            ////
////                                Credits for libraries used : GLM, glut, GLUI                        ////
////                             Big thanks to Paul Rademacher for his amazing GLUI                     ////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <GL/glut.h>
#include <glm/glm.hpp>
#include <iostream>
#include <random>
#include <stack>
using namespace glm;
using namespace std;
//PI value
const float PI = 3.1415926535f;
//The 8 vertices of the two planes
vec3 P1 = vec3(100, 100, 100);
vec3 P2 = vec3(50, -100, 100);
vec3 P3 = vec3(-50, -100, 100);
vec3 P4 = vec3(-100, 100, 100);
vec3 P5 = vec3(50, -100, -100);
vec3 P6 = vec3(100, 100, -100);
vec3 P7 = vec3(-100, 100, -100);
vec3 P8 = vec3(-50, -100, -100);

vec3 P1Normal = normalize(cross((P8 - P3), (P4 - P3))); //Normal for plane1
vec3 P2Normal = normalize(cross((P2 - P5), (P6 - P5))); //Normal for plane1

//for rotating the camera
int rotateon;
int lastx, lasty;
int xchange, ychange;
float spin = 0.0;
float spinup = 0.0;
int loop = -1;// to loop through the random guassian numbers

float t = 0; float n = 0; //Initialize current time and step number
float h = 0.001; //timestep
float g = -9.8; // Gravity

float Cr = 0.7; //Co-efficient of restitution
float Cf = 0.4; //Co-efficient of friction
vec3 vN; //Normal component of reflected velocity 
vec3 vT; //Tangential component of reflected velocity

float r = 100; //particles per second
const int NUM_PARTICLES = 5000.0; //total number of particles
float Ms = 0.06; //Mean speed of the particle
float DeltaS = 0.01; //Range of speeds of the particle


float U[100];// To store random Uniform numbers
float G[100];// To store random Guassian numbers

//Parameters for different particle generators
vec3 center = vec3(0, 0, 0); //Center point for omni directional particle generator
vec3 Vunit; //initial direction
float S0; //initial speed
vec3 VortexCenter = vec3(20.0f, 20.0f, 0.0f); //Center for a vortex
float VortexScale = 1.3; //speed for vortex
vec3 Vangle = vec3(PI, PI, PI); //angle at which vortex is rotating

//Parallel arrays for storing information of the particles
vec3 ParticlePos[NUM_PARTICLES];
vec3 ParticlePosNext[NUM_PARTICLES]; //to store position of the particle in the next time step

vec3 ParticleVel[NUM_PARTICLES];
vec3 ParticleVelNext[NUM_PARTICLES]; //to store velocity of the particles in the next time step

vec3 ParticleVelCollide[NUM_PARTICLES];
vec3 ParticleColor[NUM_PARTICLES];
float ParticleLifespan[NUM_PARTICLES];
bool ParticleActive[NUM_PARTICLES];

stack<int> inactiveStack; //Used to store the indices of inactive particles in the array
int inactiveCount; //Total number of inactive particles

//Returns returning unifromly distibuted random scalar in the range (lo....hi)
float Uniform(float lo, float hi)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(lo, hi);
    return dis(gen);
}

//Returns returning normally distributed random scalar in with the mean "mean" and standard deviation "dev"
float Guassian(float mean, float dev)
{
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<float> d(mean, dev);
    return d(gen);
}

//Returns unit vector uniformly distributed across the directions of the surface of a unit sphere.
vec3 S()
{
    float theta = Uniform(0.0, 2 * PI);
    float z = Uniform(-1.0, 1.0);
    vec3 P1;
    P1.x = sqrt(1 - pow(z, 2))*cos(theta);
    P1.y = z;
    P1.z = -sqrt(1 - pow(z, 2))*sin(theta);
    return P1;
}

//Returns unit vector uniformaly distributed about w with a maximum deflection angle theta.
vec3 Du(vec3 w, float angle)
{
    vec3 a;
    if (w.y == 0 && w.z == 0)
        a = vec3(0, 1, 0);
    else
        a = vec3(1, 0, 0);
    vec3 zaxis = w;
    vec3 xaxis = cross(a, zaxis) / length(cross(a, zaxis));
    vec3 yaxis = cross(zaxis, xaxis);
    float f = Uniform(0, 1);
    float phi = sqrt(f)*angle;
    float theta = Uniform(-PI, PI);
    vec3 Vi = vec3(cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi));
    return vec3((xaxis.x*Vi.x) + (yaxis.x*Vi.y) + (zaxis.x*Vi.z), (xaxis.y*Vi.x) + (yaxis.y*Vi.y) + (zaxis.y*Vi.z), (xaxis.z*Vi.x) + (yaxis.z*Vi.y) + (zaxis.z*Vi.z));
}

//Returns unit vector normally distributed about w with a mean "theta" and standard deviation of theta/3.
vec3 Dg(vec3 w, float angle)
{
    vec3 a;
    if (w.y == 0 && w.z == 0)
        a = vec3(0, 1, 0);
    else
        a = vec3(1, 0, 0);
    vec3 zaxis = w;
    vec3 xaxis = cross(a, zaxis) / length(cross(a, zaxis));
    vec3 yaxis = cross(zaxis, xaxis);
    float f = Guassian(0, angle / 3);
    float phi = sqrt(f)*angle;
    float theta = Uniform(-PI, PI);
    vec3 Vi = vec3(cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi));
    return vec3((xaxis.x*Vi.x) + (yaxis.x*Vi.y) + (zaxis.x*Vi.z), (xaxis.y*Vi.x) + (yaxis.y*Vi.y) + (zaxis.y*Vi.z), (xaxis.z*Vi.x) + (yaxis.z*Vi.y) + (zaxis.z*Vi.z));
}

//Returns a point uniformly distributed across the face of the disc
vec3 Cu(vec3 C, vec3 n, float R)
{
    vec3 a;
    if (n.y == 0 && n.z == 0)
        a = vec3(0, 1, 0);
    else
        a = vec3(1, 0, 0);
    vec3 zaxis = n;
    vec3 xaxis = cross(a, zaxis) / length(cross(a, zaxis));
    vec3 yaxis = cross(zaxis, xaxis);
    float f = Uniform(0, 1);
    float r = sqrt(f)*R;
    float theta = Uniform(-PI, PI);
    vec3 Vi = vec3(cos(theta)*r, sin(theta)*r, 0);
    return (vec3((xaxis.x*Vi.x) + (yaxis.x*Vi.y) + (zaxis.x*Vi.z), (xaxis.y*Vi.x) + (yaxis.y*Vi.y) + (zaxis.y*Vi.z), (xaxis.z*Vi.x) + (yaxis.z*Vi.y) + (zaxis.z*Vi.z)) + C);
}

//Returns a point normally distributed across the face of the disc
vec3 Cg(vec3 C, vec3 n, float R)
{
    vec3 a;
    if (n.y == 0 && n.z == 0)
        a = vec3(0, 1, 0);
    else
        a = vec3(1, 0, 0);
    vec3 zaxis = n;
    vec3 xaxis = cross(a, zaxis) / length(cross(a, zaxis));
    vec3 yaxis = cross(zaxis, xaxis);
    float f = Guassian(0, R / 3);
    float r = sqrt(f)*R;
    float theta = Uniform(-PI, PI);
    vec3 Vi = vec3(cos(theta)*r, sin(theta)*r, 0);
    return (vec3((xaxis.x*Vi.x) + (yaxis.x*Vi.y) + (zaxis.x*Vi.z), (xaxis.y*Vi.x) + (yaxis.y*Vi.y) + (zaxis.y*Vi.z), (xaxis.z*Vi.x) + (yaxis.z*Vi.y) + (zaxis.z*Vi.z)) + C);
}

//Efficient Random numbers : Stores arrays of random numbers from uniform and normal distribution to be used throughout the simulation 
void EfficientRandomNumbers()
{
    for (int i = 0; i < 100; i++)
    {
        U[i] = Uniform(0, 1);
        G[i] = Guassian(0, 1);
        //cout << U[i] << " , " << G[i] << endl;
    }
}

//Returns random unifrom value between Umin and Umax
float randomUniform(float Umin, float Umax, float Ui)
{
    return ((Umax - Umin)*Ui + Umin);
}

//Returns random guassian value with mean M and standard deviation D
float randomGuassian(float M, float D, float Gi)
{
    return D * Gi + M;
}

bool BarycentricTriangleOne(vec3 Point)
{
    vec3 Vn = cross((P8 - P3), (P7 - P8));
    float TwoA = length(Vn);
    vec3 n = Vn / TwoA;
    float u = dot(cross((P7 - P8), (Point - P8)), n / TwoA);
    float v = dot(cross((P3 - P7), (Point - P7)), n / TwoA);
    float w = 1 - u - v;

    if (u >= 0 && v >= 0 && u + v <= 1) { return(true); }

}

bool BarycentricTriangleTwo(vec3 Point)
{
    vec3 Vn = cross((P4 - P3), (P7 - P4));
    float TwoA = length(Vn);
    vec3 n = Vn / TwoA;
    float u = dot(cross((P7 - P4), (Point - P4)), n / TwoA);
    float v = dot(cross((P3 - P7), (Point - P7)), n / TwoA);
    float w = 1 - u - v;

    if (u >= 0 && v >= 0 && u + v <= 1) { return(true); }
}

bool InsidePolygonCheck(vec3 Point)
{
    bool InsideOne = BarycentricTriangleOne(Point);
    bool InsideTwo = BarycentricTriangleTwo(Point);
    if (InsideOne == TRUE || InsideTwo == TRUE)
    {
        return(true);
    }
    else
    {
        return(false);
    }
}

//Assigning paremeters to the particles for different particle generators
vec3 Directed = normalize(P3 - center); //To store random Direction for our Directed particle generator
void ParticleParameters(int Ptype)
{

    for (int i = 0; i < NUM_PARTICLES; i++)
    {

        switch (Ptype)
        {
        default:
            ParticlePos[i] = center;
            ParticlePosNext[i] = center;
            Vunit = S();
            loop = loop + 1;
            if (loop > 99)loop = 0;
            S0 = randomGuassian(Ms, DeltaS / 3, G[loop]);
            ParticleVel[i] = S0 * Vunit;
            ParticleVelNext[i] = S0 * Vunit;
            ParticleLifespan[i] = abs(randomGuassian(3.0, 1, G[loop]));
            if (ParticleLifespan[i] <= 1.0f) { ParticleLifespan[i] = 1; }
            else if (ParticleLifespan[i] >= 4.0f) { ParticleLifespan[i] = 4; }
            ParticleColor[i] = vec3(U[rand() % 100], U[rand() % 100], U[rand() % 100]);

            ParticleActive[i] = true;
            //cout << "Speed" << S0 << " | Velocity:" << ParticleVel[i].x << " , " << ParticleVel[i].y << " , " << ParticleVel[i].z << endl;
            //cout << "Speed" << S0 << " | Velocity:" << ParticleVelNext[i].x << " , " << ParticleVelNext[i].y << " , " << ParticleVelNext[i].z << endl;
            //cout << i << " Lifespan: " << ParticleLifespan[i] << " | Gloop: " << G[loop] << " Loop: " << loop << endl;
            //cout << "Omnidirectional particles" << endl;
            break;

        case 1:
            ParticlePos[i] = center;
            Vunit = Du(Directed, PI / 8);
            S0 = randomGuassian(Ms, DeltaS / 3, G[loop]);
            loop = loop + 1;
            if (loop > 99)loop = 0;
            ParticleVel[i] = S0 * Vunit;
            ParticleLifespan[i] = abs(randomGuassian(3.0, 1, G[loop]));
            if (ParticleLifespan[i] <= 1.0f) { ParticleLifespan[i] = 1; }
            else if (ParticleLifespan[i] >= 4.0f) { ParticleLifespan[i] = 4; }
            ParticleColor[i] = vec3(U[rand() % 100], U[rand() % 100], U[rand() % 100]);
            ParticleActive[i] = true;
            //cout << "Directed particles" << endl;
            break;
        case 2:
            ParticlePos[i] = center;
            Vunit = Du(Directed, PI / 3);
            S0 = randomGuassian(Ms, DeltaS / 3, G[loop]);
            loop = loop + 1;
            if (loop > 99)loop = 0;
            ParticleVel[i] = S0 * Vunit;
            ParticleLifespan[i] = abs((3.0, 0.4, G[loop]));
            ParticleActive[i] = true;
            cout << "Directed particles" << endl;
            break;
        }
    }
}

void GenerateParticles(int N, int Ptype)
{
    switch (Ptype)
    {
    default:
        ParticlePos[N] = center;
        ParticlePosNext[N] = center;
        Vunit = S();
        loop = loop + 1;
        if (loop > 99)loop = 0;
        S0 = randomGuassian(Ms, DeltaS / 3, G[loop]);
        ParticleVel[N] = S0 * Vunit;
        ParticleVelNext[N] = S0 * Vunit;
        ParticleLifespan[N] = abs(randomGuassian(3.0, 1, G[loop]));
        if (ParticleLifespan[N] <= 1.0f) { ParticleLifespan[N] = 1; }
        else if (ParticleLifespan[N] >= 4.0f) { ParticleLifespan[N] = 4; }
        ParticleColor[N] = vec3(U[loop], U[loop], U[loop]);

        ParticleActive[N] = true;
        inactiveStack.pop();
        break;
    case 1:
        ParticlePos[N] = center;
        ParticlePosNext[N] = center;
        Vunit = Du(Directed, PI / 3);
        S0 = randomGuassian(Ms, DeltaS / 3, G[loop]);
        loop = loop + 1;
        if (loop > 99)loop = 0;
        ParticleVel[N] = S0 * Vunit;
        ParticleVelNext[N] = S0 * Vunit;
        ParticleLifespan[N] = abs(randomGuassian(3.0, 1, G[loop]));
        if (ParticleLifespan[N] <= 1.0f) { ParticleLifespan[N] = 1; }
        else if (ParticleLifespan[N] >= 4.0f) { ParticleLifespan[N] = 4; }
        ParticleColor[N] = vec3(U[loop], U[loop], U[loop]);
        ParticleActive[N] = true;

        inactiveStack.pop();
        break;
    }
}

void testAndDeactivate()
{
    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        if (ParticleActive[i] == true)
        {

            if (abs(length(center - ParticlePos[i])) > 200)
            {
                inactiveStack.push(i);
                inactiveCount += 1;
                ParticleActive[i] = false;
                GenerateParticles(inactiveStack.top(), 0);
            }
        }
    }
}

vec3 VortexVelocity(vec3 Point)
{
    vec3 D = vec3(Point - VortexCenter);
    vec3 V = cross(Vangle, D);
    return(V);
}

float factor(vec3 Vector)
{
    vec3 D = vec3(Vector - VortexCenter);
    vec3 V = cross(Vangle, D);
    float f1 = 1 + pow(D.x, 2) + pow(D.y, 2) + pow(D.z, 2);
    return 1 / (f1 + VortexScale);
}

void CollisionCheck(float H) //Checking if the point has hit the left infinite plane
{
    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        if (ParticleActive[i] == true)
        {
            ParticleVelNext[i].x = ParticleVel[i].x;
            ParticleVelNext[i].y = ParticleVel[i].y + (g * H);
            ParticleVelNext[i].z = ParticleVel[i].z;
            //Adding Vortex effect
            ParticleVelNext[i].x += ((VortexVelocity(ParticlePos[i]).x) - ParticleVelNext[i].x)*factor(ParticlePos[i]);
            ParticleVelNext[i].y += ((VortexVelocity(ParticlePos[i]).y) - ParticleVelNext[i].y)*factor(ParticlePos[i]);
            ParticleVelNext[i].z += ((VortexVelocity(ParticlePos[i]).z) - ParticleVelNext[i].z)*factor(ParticlePos[i]);

            ParticlePosNext[i] = ParticlePos[i] + ParticleVelNext[i];
            float dN = dot((ParticlePos[i] - P3), P1Normal);
            float dN1 = dot((ParticlePosNext[i] - P3), P1Normal);
            if (dN*dN1 < 0)
            {
                ParticlePosNext[i] = ParticlePosNext[i] - ((1.0f + Cr)*dN1*P1Normal);
                if (InsidePolygonCheck(ParticlePosNext[i]) == true) //Checking if the point lies wihtin the finite plane
                {
                    vN = dot((ParticleVelNext[i]), P1Normal)*P1Normal;
                    vT = ParticleVelNext[i] - vN;
                    ParticleVelCollide[i] = ((-Cr)*vN + (1 - Cf)*vT);
                    ParticleVelNext[i] = ParticleVelCollide[i];
                }

            }
            ParticlePos[i] = ParticlePosNext[i];
        }
    }

}



void integrate(float H)
{
    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        if (ParticleActive[i] == true)
        {
            //ParticlePos[i] = ParticlePosNext[i];
            //ParticleVel[i] = ParticleVelNext[i];
            ParticleVelNext[i].x = ParticleVel[i].x;
            ParticleVelNext[i].y = ParticleVel[i].y + (g * H);
            ParticleVelNext[i].z = ParticleVel[i].z;
            ParticlePos[i] = ParticlePos[i] + ParticleVelNext[i];

        }

    }
}

void init(void)
{
    //glClearColor(0, 0, 0, 0.0);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glShadeModel(GL_SMOOTH);
    //glEnable(GL_POINT_SMOOTH);
    //glEnable(GL_LIGHTING);
    //glEnable(GL_LIGHT0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, 1.0, 1, 600);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();


    glPointSize(0.2);
    // Set eye point and lookat point
    gluLookAt(0, 200, 300, 0, 0, 0, 0, 1, 0);
}


void reshapeFunc(GLint newWidth, GLint newHeight)
{
    if (newWidth > newHeight) // Keep a square viewport
        glViewport((newWidth - newHeight) / 2, 0, newHeight, newHeight);
    else
        glViewport(0, (newHeight - newWidth) / 2, newWidth, newWidth);
    init();
    glutPostRedisplay();
}


void mouse(int button, int state, int x, int y)
{
    switch (button) {
    case GLUT_LEFT_BUTTON:
        if (state == GLUT_DOWN) {
            lastx = x;
            lasty = y;
            xchange = 0;
            ychange = 0;
            rotateon = 1;
        }
        else if (state == GLUT_UP) {
            xchange = 0;
            ychange = 0;
            rotateon = 0;
        }
        break;

    default:
        break;
    }
}

void motion(int x, int y)
{
    xchange = x - lastx;
    ychange = y - lasty;
}


void idleCallback(void)
{
    if (t < 80)
    {
        testAndDeactivate();
        CollisionCheck(h);
        //integrate(h);


        n = n + 1;
        t = n * h;
    }




    if (rotateon) {
        spin += xchange / 250.0;
        if (spin >= 360.0) spin -= 360.0;
        if (spin < 0.0) spin += 360.0;
        spinup -= ychange / 250.0;
        if (spinup > 89.0) spinup = 89.0;
        if (spinup < -89.0) spinup = -89.0;
    }
    glutPostRedisplay();
}

void display(void)
{
    GLfloat lines[] = { 1, 1, 1 };

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glPushMatrix();

    //rotate the view
    glRotatef(spinup, 1.0, 0.0, 0.0);
    glRotatef(spin, 0.0, 1.0, 0.0);


    glBegin(GL_LINES);
    glMaterialfv(GL_FRONT, GL_AMBIENT, lines);
    glColor3f(1.0, 1.0, 1.0);
    glVertex3f(P1.x, P1.y, P1.z); glVertex3f(P2.x, P2.y, P2.z);
    glVertex3f(P3.x, P3.y, P3.z); glVertex3f(P4.x, P4.y, P4.z);
    glVertex3f(P8.x, P8.y, P8.z); glVertex3f(P7.x, P7.y, P7.z);
    glVertex3f(P6.x, P6.y, P6.z); glVertex3f(P5.x, P5.y, P5.z);
    glVertex3f(P1.x, P1.y, P1.z); glVertex3f(P6.x, P6.y, P6.z);
    glVertex3f(P2.x, P2.y, P2.z); glVertex3f(P5.x, P5.y, P5.z);
    glVertex3f(P4.x, P4.y, P4.z); glVertex3f(P7.x, P7.y, P7.z);
    glVertex3f(P3.x, P3.y, P3.z); glVertex3f(P8.x, P8.y, P8.z);




    glEnd();
    glPopMatrix();

    glPushMatrix();
    glRotatef(spinup, 1.0, 0.0, 0.0);
    glRotatef(spin, 0.0, 1.0, 0.0);
    glBegin(GL_POINTS);
    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        glColor3f(ParticleColor[i].x, ParticleColor[i].y, ParticleColor[i].z);
        glVertex3f(ParticlePos[i].x, ParticlePos[i].y, ParticlePos[i].z);

    }
    glEnd();
    /*for (int i = 0; i < NUM_PARTICLES; i++)
    {
        glTranslatef(ParticlePos[i].x, ParticlePos[i].y, ParticlePos[i].z);
        glutSolidSphere(1.2, 10, 10);
    }*/

    glPopMatrix();
    glutSwapBuffers();
}

int main(int argc, char** argv)
{

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(700, 700);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Particle System");
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glClearColor(0.0, 0.0, 0.0, 0.0);
    init();
    rotateon = 0;
    EfficientRandomNumbers();
    ParticleParameters(0);

    glutDisplayFunc(display);
    glutIdleFunc(idleCallback);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutReshapeFunc(reshapeFunc);
    glutMainLoop();
    return 0;
}
