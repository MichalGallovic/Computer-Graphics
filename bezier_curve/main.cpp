// Michal Gallovic
// Computer Graphics
// Draw bezier curve, Move the ball over the surface of the curve

#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <vector>
#include <math.h>
#include <sstream>
#include <string>
#include <stdio.h>
#include <iostream>
#if defined(__APPLE__)
#include "glm/glm.hpp"
#else
#include <glm/glm.hpp>
#endif

using namespace glm;
using namespace std;

// HELPERS
string convertInt(int number) {
    stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
}
void renderBitmapString(float x, float y, void *font, string s) {
    
    glRasterPos2f(x,y);
    for (int i=0; i < s.size(); i++)
    {
        glutBitmapCharacter(font, s[i]);
    }
}
//#######
float krok = 0.01;


class Ball{
public:
    
    Ball(){
        s = 0;
        r = 0.1;
        v = vec3(0,0,0);
        a = vec3(0,-9.81,0);
    }
    
    void state_change(float deltat){
        
        
        C += v*deltat;

        v += a*deltat;
        
        s += glm::length(v*deltat);
        
    }

    void draw(){
        glColor3f(1,1,1);
        glTranslatef(C.x,C.y,C.z);
        
        if(this->isGoingToRight){
            glRotatef(-s/6.28/r*360,0,0,1);
            glutSolidSphere(r,10,10);
            glRotatef(s/6.28/r*360,0,0,1);
        }else{
            glRotatef(s/6.28/r*360,0,0,1);
            glutSolidSphere(r,10,10);
            glRotatef(-s/6.28/r*360,0,0,1);
            
        }
        glTranslatef(-C.x,-C.y,-C.z);

    }

    vec3 C,v,a;
    float r,s;
    bool isGoingToRight;
};


class Curve{
public:
    Curve(){
        active_P = -1;
//        this->ball = ball;
//        points.push_back(vec3(-0.8,0.8,0.0));
//        points.push_back(vec3(-0.6,-0.8,0.0));
//        points.push_back(vec3(0.6,-0.8,0.0));
//        points.push_back(vec3(0.8,0.8,0.0));
//        this->generateCenterPoints();
    }
    void draw(){
        draw_curve();
        draw_points();
//        draw_center_points();
    }
    void add_point(vec3 newP){
        points.push_back(newP);
        this->generateCenterPoints();
    }
    
    bool isZeroVector (vec3 p) {
        return p.x == 0 && p.y == 0 && p.z == 0;
    }
    
    void generateCenterPoints() {
        this->bezierCenterPoints.clear();
        vec3 centerPoint;
        vec3 P1,P2, tangent, crossProduct;
        
        for (float t = 0; t <= 1.0; t+=0.001) {
            P2 = this->Bezier(t+0.001);
            P1 = this->Bezier(t-0.001);
            
            centerPoint = this->Bezier(t);
            tangent = glm::normalize(P1 - P2);
            crossProduct = cross(tangent, vec3(0,0,1));
            
            if (!glm::all(glm::isnan(crossProduct)) && this->ball != NULL) {
                centerPoint += normalize(crossProduct) * this->ball->r;
            }
            
            bezierCenterPoints.push_back(centerPoint);
        }
        
    }
    void picking(vec3 p){
        float dist_min=1000, index  = -1;
        for(int i = 0; i < points.size(); i++)
        {
            if( glm::distance(p,points[i]) < dist_min)
            {
                index = i;
                dist_min = glm::distance(p,points[i]);
            }
        }
        if(dist_min < 0.02)active_P = index;
        else active_P = -1;
    }
    void point_move(vec3 p){
        if(active_P > -1 && active_P < points.size()) {
            points[active_P] = p;
            this->generateCenterPoints();
        }
    }
    
    void draw_points(){
        glPointSize(5);
        glBegin(GL_POINTS);
        for(int i=0; i < points.size(); i++)
        {
            if(active_P == i)glColor3f(1,0,0);
            else glColor3f(0,1,0);
            glVertex3f(points[i].x, points[i].y, points[i].z);
        }
        glEnd();
    }
    
    void draw_center_points() {
        glPointSize(2);
        glColor3f(0, 1, 0);
        glBegin(GL_POINTS);
        for(int i=0; i < bezierCenterPoints.size(); i++)
        {
            glVertex3f(bezierCenterPoints[i].x, bezierCenterPoints[i].y, bezierCenterPoints[i].z);
        }
        glEnd();
    }
    
    vec3 Casteljau(int i, int j,float t){
        if(i == 0) return points[j];
        return (1-t)*Casteljau(i-1,j-1,t)+t*Casteljau(i-1,j,t);
    }
    vec3 Bezier(float t){
        int n = (int)points.size() - 1;
        if(n < 1) return vec3(0,0,0);
        return Casteljau(n,n,t);
    }
    
    float clossest(vec3 mouse){
        float dist_min = 1, clossest = -1;
        for (float t=0; t<1.0;t+=krok){
            float distance = glm::distance(mouse, Bezier(t));
            if(distance<dist_min){
                dist_min= distance;
                clossest = t;
            }
        }
        return clossest;
    }
    
    void draw_curve(){
        vec3 point;
        glPointSize(5);
        
        //draw Bezier
        glColor3f(0, 0, 1);
        glBegin(GL_LINE_STRIP);
        for (float t = 0.0; t<=1; t+=krok){
            
            point = Bezier(t);
            glVertex3f(point.x,point.y,point.z);
        }
        
        glEnd();
        
        glColor3f(1,0,0);

        
    }
    
    void drawTangent(vec3 mousePos) {
        if (points.size() > 3) {
            float t = clossest(mousePos);
            vec3 point1 =  Casteljau(points.size()-2, points.size()-2, t);
            vec3 point2 = Casteljau(points.size()-2, points.size()-1, t);
            glColor3f(1, 0, 0);
            glBegin(GL_LINES);
            glVertex3f(point1.x, point1.y, point1.z);
            glVertex3f(point2.x, point2.y, point2.z);
            glEnd();
        }
    }
    
    vector<vec3> points;
    vector<vec3> bezierCenterPoints;
    Ball *ball;
    int active_P;
    
};

class Camera {
public:
    void setOrthographicProjection(int w, int h) {
        glMatrixMode(GL_PROJECTION);
        glViewport(0, 0, w, h);
        glPushMatrix();
        glLoadIdentity();
        gluOrtho2D(0, w, 0, h);
        glMatrixMode(GL_MODELVIEW);
    }
    void resetPerspectiveProjection() {
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glLoadIdentity();
        glMatrixMode(GL_MODELVIEW);
    }
};

class Scene {
public:
    Scene(){
        //basic scene settings
        glClearColor(0.0, 0.0, 0.0, 1.0);
        menuIsVisible = true;
        startSimulation = false;
        timebase = 0;
        
        camera = new Camera();
        curve = new Curve();
        ball = new Ball();

    }
    
    void iter(){
        float t;
        
        vec3 P1,P2,tangent, clossestBezierPoint;
        int clossetsPointIndex = 0;
        ball->state_change(0.01);
        vec3 clossestCenterPoint = vec3(100,100,100);
        t = curve->clossest((vec3)ball->C);
        P2 = curve->Bezier(t+0.001f);
        P1 = curve->Bezier(t-0.001f);
        clossestBezierPoint = curve->Bezier(t);
        tangent = glm::normalize(P2-P1);
        
//        for (int i = 0; i < curve->bezierCenterPoints.size(); i++) {
//            if (distance(curve->bezierCenterPoints[i], ball->C) < distance(clossestCenterPoint,ball->C)) {
//                clossestCenterPoint = curve->bezierCenterPoints[i];
//                clossetsPointIndex = i;
//            }
//        }
//        P1 = curve->bezierCenterPoints[clossetsPointIndex];
//        P2 = curve->bezierCenterPoints[clossetsPointIndex + 1];
//        tangent = glm::normalize(P2-P1);
        
        if(glm::dot(tangent,ball->v) > 0.0){

            ball->isGoingToRight = true;
            
            float newChange = (tangent*glm::length(ball->v) - 0.02f).x;
//            printf("%.32f\n",(abs(lastChange) - abs(newChange)));
            if((abs(lastChange) - abs(newChange)) != 0.0) {
//                printf("vbefore: %f %f\n", ball->v.x, ball->v.y);
                if (this->isInLocalMin) {
//                    ball->v = tangent*glm::length(ball->v)- 1.f;
                    if ((totalframes - localMinFrame)/60 > 0.5) {
                        this->resetIter();
                        this->spawnBall();
                    }

                } else {
                    ball->v = tangent*glm::length(ball->v)- 0.02f;
                }
                
//                printf("tangentChange: %f %f\n",(-tangent*glm::length(ball->v) + 0.01f).x,(-tangent*glm::length(ball->v) + 0.01f).y);
                ball->isGoingToRight = true;
//                printf("vafter: %f %f\n", ball->v.x, ball->v.y);
                lastChange = ball->v.x;
            } else {
//                ball->v = -tangent*glm::length(ball->v/100.f);
                this->isInLocalMin = true;
                this->localMinFrame = totalframes;
                printf("up");
            }

        }else{

            float newChange = (-tangent*glm::length(ball->v) + 0.02f).x;
//            printf("%.32f\n",(abs(lastChange) - abs(newChange)));
            if((abs(lastChange) - abs(newChange)) != 0.0) {
//                printf("vbefore: %f %f\n", ball->v.x, ball->v.y);
            
                if (this->isInLocalMin) {
//                    ball->v = -tangent*glm::length(ball->v)+ 1.f;
                    if ((totalframes - localMinFrame)/60 > 0.5) {
                        this->resetIter();
                        this->spawnBall();
                    }

                } else {
                    ball->v = -tangent*glm::length(ball->v)+ 0.02f;
                }
//                printf("tangentChange: %f %f\n",(-tangent*glm::length(ball->v) + 0.01f).x,(-tangent*glm::length(ball->v) + 0.01f).y);
                ball->isGoingToRight = false;
//                printf("vafter: %f %f\n", ball->v.x, ball->v.y);
                lastChange = ball->v.x;
            } else {
//                ball->v = tangent*glm::length(ball->v/100.f);
                this->isInLocalMin = true;
                this->localMinFrame = totalframes;
                printf("down");
            }
        }
        
        clossestCenterPoint = dvec3(100,100,100);
        
        for (int i = 0; i < curve->bezierCenterPoints.size(); i++) {
            if (distance(curve->bezierCenterPoints[i], ball->C) < distance(clossestCenterPoint,ball->C)) {
                clossestCenterPoint = curve->bezierCenterPoints[i];
            }
        }


        if (clossestCenterPoint != curve->bezierCenterPoints[0] && clossestCenterPoint != curve->bezierCenterPoints[curve->bezierCenterPoints.size()-1]) {
            ball->C = clossestCenterPoint;
        }
        
        
    }
    
    void draw() {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glLoadIdentity();
        this->drawText();
        curve->draw();
        if (this->startSimulation) {
            ball->draw();
        }
        glutSwapBuffers();
    }
    
    void drawText() {
        glColor3d(1.0, 0.0, 0.0);
        camera->setOrthographicProjection(w,h);
        glPushMatrix();
        glLoadIdentity();
        
        if (menuIsVisible) {
            drawMenu();
        } else {
            drawFPS();
            if(totalframes / 60 == 3) {
                hintsVisible = false;
            }
            if (hintsVisible) {
                drawHints();
            } else {
                minimizeHints();
            }
        }
        
        glPopMatrix();
        camera->resetPerspectiveProjection();
    }
    
    void drawMenu() {
        renderBitmapString(w/2 - 50, h/2, GLUT_BITMAP_HELVETICA_18, "Bezier Curve");
        renderBitmapString(w/2 - 80, h/2 - 30, GLUT_BITMAP_HELVETICA_18, "Press Enter to Start");
        
    }
    void drawHints() {
        renderBitmapString(10, h - 20, GLUT_BITMAP_HELVETICA_12, "Draw curve with mouse clicks");
        renderBitmapString(10, h - 40, GLUT_BITMAP_HELVETICA_12, "Add point - Right mouse button");
        renderBitmapString(10, h - 60, GLUT_BITMAP_HELVETICA_12, "Move point - Left mouse button");
        renderBitmapString(10, h - 80, GLUT_BITMAP_HELVETICA_12, "Press s to release ball");
        renderBitmapString(10, h - 100, GLUT_BITMAP_HELVETICA_12, "Press ESC to quit");
    }
    void minimizeHints() {
        renderBitmapString(10, h - 20, GLUT_BITMAP_HELVETICA_12, "For Hints press & hold h");
    }
    void drawFPS() {
        string s_fps = "FPS: " + convertInt(fps);
        renderBitmapString(w - 60, h - 20, GLUT_BITMAP_HELVETICA_12, s_fps);
    }
    void countFPS() {
        frame++;
        totalframes++;
        time=glutGet(GLUT_ELAPSED_TIME);
        
        if (time - timebase > 1000) {
            fps = frame*1000.0/(time-timebase);
            timebase = time;
            frame = 0;
        }
    }

    void reshape(int width, int height) {
        w = width;
        h = height;
        camera->setOrthographicProjection(width, height);
    }
    
    void mouseClick(int button, int state,int x, int y){
        if(button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
        {
            curve->add_point(image_to_coor(x,y));
        }
        if(button == GLUT_RIGHT_BUTTON && state == GLUT_UP)
        {
            curve->active_P = -1;
        }
        
        if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
        {
            curve->picking(image_to_coor(x,y));
        }
        if(button == GLUT_LEFT_BUTTON && state == GLUT_UP)
        {
            curve->active_P = -1;
        }
    }
    void mouseDrag(int x, int y){
        curve->point_move(image_to_coor(x,y));
    }
    vec3 image_to_coor(int x, int y){
        int width,height;
        float xre,yre;
        
        width = glutGet(GLUT_WINDOW_WIDTH);
        height = glutGet(GLUT_WINDOW_HEIGHT);
        
        xre = 2*(float)x/width - 1;
        yre = -2*(float)y/height + 1;

        return vec3(xre,yre,0);
    }

    void spawnBall() {
        ball->C = curve->bezierCenterPoints[0];
    }
    
    void resetIter() {
        ball = new Ball();
        this->isInLocalMin = false;
    }
    
    bool hintsVisible, menuIsVisible, startSimulation, isInLocalMin;
    int localMinFrame;
    int w,h;
    int frame, time, timebase;
    int totalframes;
    float fps;
    float lastChange;
    Camera *camera;
    Curve *curve;
    Ball *ball;
};




// GLOB VARS
Scene *scene;

//

void display(void) {
    scene->countFPS();
    scene->draw();
    if (scene->startSimulation) {
        scene->iter();
    }
}

void repeat() {
    
    glutPostRedisplay();
}

void reshape(int width, int height) {
    scene->reshape(width, height);
}

void keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case 'h':
            scene->hintsVisible = true;
            break;
        case 13:
            scene->menuIsVisible = false;
            scene->hintsVisible = true;
            //reset scene total frames
            scene->totalframes = 0;
            break;
        case 27:
            exit(0);
            break;
        case 's':
            scene->startSimulation = true;
            scene->resetIter();
            scene->spawnBall();
            break;
    }
}

void keyboardUp(unsigned char key, int x, int y) {
    switch (key) {
        case 'h':
            scene->hintsVisible = false;
            break;
    }
}

void mouse(int button, int state, int x, int y)
{
    scene->mouseClick(button,state,x,y);
}
void motion(int x, int y)
{
    scene->mouseDrag(x,y);
}

void init() {
    scene = new Scene();
}

int main(int argc, char *argv[])
{
    glutInit(&argc, argv);
    glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowSize(600, 600);
    glutCreateWindow("Bezier Curves");
    init();
    glutDisplayFunc(display);
    glutIdleFunc(repeat);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutKeyboardUpFunc(keyboardUp);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutMainLoop();
    
    return 0;
}
