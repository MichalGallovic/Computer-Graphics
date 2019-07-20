#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <iostream>
#include <vector>
#include <map>
#include <string>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <sstream>
#include <stdio.h>

using namespace std;
using namespace glm;


const vec3 nullvec = vec3(0,0,0);
const vec3 red = vec3(1,0,0);
const vec3 green = vec3(0,1,0);
const vec3 blue = vec3(0,0,1);
const vec3 yellow = vec3(1,1,0);
const vec3 magenta = vec3(1,0,1);
const vec3 azure = vec3(0,1,1);
const vec3 brown = vec3((float) 160/255,(float) 82/255,(float) 45/255);
const vec3 gray = vec3((float) 240/255,(float) 250/255,(float) 240/255);
const vec3 black = vec3(0,0,0);
const vec3 white = vec3(1,1,1);
const GLuint KOCH_SNOWFLAKE = 1;
const GLuint SQUARE_CARPET = 2;
const GLuint SIERPINSKY_TRIANGLE = 3;
const GLuint HILBERT_CURVE = 4;


//helpers
vector<vec3>hilbertPoints;

vec2 normalizePosition(int x, int y) {
    float width = (float)glutGet(GLUT_WINDOW_WIDTH);
    float height = (float)glutGet(GLUT_WINDOW_HEIGHT);
    
    float pozX =  (x/(width/2))-1;
    float pozY =  (-1)*(y/(height/2))+1;
    vec2 point = vec2(pozX, pozY);
    
    return point;
}

class Ball {
public:
    Ball(vec3 center = nullvec, float ballSize = 0.1, vec3 color = blue) {
        this->pos = center;
        this->ballSize = ballSize;
        this->color = color;
        
        this->outerBox = 1.0;
        this->innerBox = outerBox - ballSize;
    }
    
    void draw() {
        vec4 alphaColor = vec4(blue, 0.0f);
        glColor4fv((GLfloat *) &alphaColor);
        glTranslatef(pos.x,pos.y,pos.z);
        glutSolidSphere(ballSize,20,20);
    }
    
    void move() {
        if((pos.x+vel.x) >= innerBox){
            pos.x = pos.x - 2*(pos.x+vel.x - innerBox);
            vel.x = -vel.x;
        } else if((pos.x-vel.x) <= -innerBox) {
            pos.x = pos.x - 2*(pos.x+vel.x - (-innerBox));
            vel.x = -vel.x;
        }
        pos.x+=vel.x;
        
        if((pos.y+vel.y) >= innerBox){
            pos.y = pos.y - 2*(pos.y+vel.y - innerBox);
            vel.y = -vel.y;
        } else if((pos.y-vel.y) <= -innerBox) {
            pos.y = pos.y - 2*(pos.y+vel.y - (-innerBox));
            vel.y = -vel.y;
        }
        pos.y+=vel.y;
        
        if((pos.z+vel.z) >= innerBox){
            pos.z = pos.z - 2*(pos.z+vel.z - innerBox);
            vel.z = -vel.z;
        } else if((pos.z-vel.z) <= -innerBox) {
            pos.z = pos.z - 2*(pos.z+vel.z - (-innerBox));
            vel.z = -vel.z;
        }
        pos.z+=vel.z;
    }
    
    void setColor(vec3 color) {
        this->color = color;
    }
    
    void setPos(vec3 pos){
        this->pos = pos;
    }
    vec3 getPos(void){
        return this->pos;
    }
    
    void setVel(vec3 vel){
        this->vel = vel;
    }
    
    vec3 getVel(void){
        return this->vel;
    }
    
    void setSize(float ballSize){
        this->ballSize = ballSize;
        this->innerBox = innerBox - ballSize;
    }
    
    float getSize(void){
        return this->ballSize;
    }
    
    void setOuterBox(float box){
        this->outerBox = box;
        this->innerBox = outerBox - ballSize;
    }
private:
    vec3 pos, vel, color;
    float ballSize;
    float outerBox;
    float innerBox;
    
};

class Fractal {

public:
    static void SierpinskyLS(vec3 point,float size, int level){
        if(level == 1)
        {
            triangle(point,size);
        }
        else
        {
            SierpinskyLS(point,size/2,level-1);
            SierpinskyLS(point+vec3(size/2,0,0),size/2,level-1);
            SierpinskyLS(point+vec3(0,size/2,0),size/2,level-1);
        }
    }
    
    static void SierpinskySquareCarpet(vec3 point, float size, int level) {
        if(level != 0) {
            square(point+vec3(size/3,size/3,0), size/3);
            SierpinskySquareCarpet(point, size/3, level -1);
            SierpinskySquareCarpet(point+vec3(size/3,0,0), size/3, level -1);
            SierpinskySquareCarpet(point+vec3(2*(size/3),0,0), size/3, level -1);
            SierpinskySquareCarpet(point+vec3(0,size/3,0), size/3, level -1);
            SierpinskySquareCarpet(point+vec3(2*(size/3),size/3,0), size/3, level -1);
            SierpinskySquareCarpet(point+vec3(0,2*(size/3),0), size/3, level -1);
            SierpinskySquareCarpet(point+vec3(size/3,2*(size/3),0), size/3, level -1);
            SierpinskySquareCarpet(point+vec3(2*(size/3),2*(size/3),0), size/3, level -1);
        }
    }
    //not fractal but...seems nice here
    static void chessFloor(vec3 point, double size, int count) {
        int colorToggle = false;
        vec4 alphaWhite = vec4(white,0.0f);
        vec4 alphaBlack = vec4(black, 0.0f);

        // for precision purposes added round - when /7 it made aditional lines
        for (double y = point.y; round(y*1000)/1000.0 < point.y + size; y += size/count) {
            
            for (double x = point.x; round(x*1000)/1000.0 < point.x + size; x += size/count) {

                if (colorToggle) {
                    glColor4fv((GLfloat *)&alphaWhite);
                } else {
                    glColor4fv((GLfloat *)&alphaBlack);
                }
                square(vec3(x,y,point.z), size/count);
                colorToggle = !colorToggle;
            }
            
            if (count % 2 == 0) {
                colorToggle = !colorToggle;
            }
        }
    }
    
    static void SquareCarpet(vec3 point, double size, int level) {
        if(level !=0) {
            square(point + vec3(size/3,size/3,0), size/3);
            SquareCarpet(point, size/2, level-1);
            SquareCarpet(point+vec3(size/2,0,0), size/2, level-1);
            SquareCarpet(point+vec3(0,size/2,0), size/2, level-1);
            SquareCarpet(point+vec3(size/2,size/2,0), size/2, level-1);
        }
    }

    static void HilbertCurve(float x, float y,float z,float xi, float xj, float yi, float yj, int n) {
        if (n <= 0) {
            hilbertPoints.push_back(vec3(x + (xi + yi)/2, y + (xj + yj)/2,z));

        } else {
            HilbertCurve(x,           y, z,           yi/2, yj/2,  xi/2,  xj/2, n-1);
            HilbertCurve(x+xi/2,      y+xj/2 , z,     xi/2, xj/2,  yi/2,  yj/2, n-1);
            HilbertCurve(x+xi/2+yi/2, y+xj/2+yj/2, z, xi/2, xj/2,  yi/2,  yj/2, n-1);
            HilbertCurve(x+xi/2+yi,   y+xj/2+yj, z,  -yi/2,-yj/2, -xi/2, -xj/2, n-1);
        }
    }
    
    static void KochSnowFlake(int level,float size){
        
        if(level == 1)
        {
            glBegin(GL_LINES);
            glVertex3f(0,0,0);
            glVertex3f(3*size,0,0);
            glEnd();
        }
        else
        {
            glPushMatrix();
            KochSnowFlake(level-1,size/3);
            glTranslatef(size,0,0);
            glRotatef(60,0,0,1);
            KochSnowFlake(level-1,size/3);
            glTranslatef(size,0,0);
            glRotatef(-120,0,0,1);
            KochSnowFlake(level-1,size/3);
            glTranslatef(size,0,0);
            glRotatef(60,0,0,1);
            KochSnowFlake(level-1,size/3);
            glPopMatrix();
        }
    }
    
private:
    static void triangle(vec3 a, float size){
        vec3 b = a + vec3(size,0,0), c = a + vec3(0,size,0);
        glBegin(GL_POLYGON);
        glVertex3fv((GLfloat*)&a);
        glVertex3fv((GLfloat*)&b);
        glVertex3fv((GLfloat*)&c);
        glEnd();
    }
    
    static void square(vec3 a, float size) {
        vec3 b = a + vec3(size,0,0);
        vec3 c = b + vec3(0,size,0);
        vec3 d = c - vec3(size,0,0);
        glBegin(GL_POLYGON);
        glVertex3fv((GLfloat*)&a);
        glVertex3fv((GLfloat*)&b);
        glVertex3fv((GLfloat*)&c);
        glVertex3fv((GLfloat*)&d);
        glEnd();
    }
    
    
};

class Room {
public:
    Room(vec3 center = nullvec, float roomSize = 1.0) {
        this->center = center;
        this->roomSize = roomSize;
        this->wallColors["top"] = black;
        this->wallColors["bottom"] = brown;
        this->wallColors["left"] = gray;
        this->wallColors["right"] = gray;
        this->wallColors["far"] = gray;
        this->wallColors["near"] = gray;
        this->timer = clock();
        
        this->chessBoardTiles = 5;
        this->fractalsCompiled = false;
        
        
    }
    
    void setRoomSize(float size) {
        this->roomSize = size;
    }
    
    void incrBoardCount() {
        chessBoardTiles++;
    }
    
    void decrBoardcount() {
        if(chessBoardTiles != 0) {
            chessBoardTiles--;
        }
    }
    
    void draw() {
        cube(center, roomSize);
        
    }
    
    float getRoomSize() {
        return roomSize;
    }
    
    vec3 getCenter() {
        return center;
    }
    
private:
    
    vec3 center;
    float roomSize;
    map<string, vec3> wallColors;
    clock_t timer;
    int chessBoardTiles;
    bool fractalsCompiled;
    
    void cubeSide(vec3 center, float size, vec3 color){
        glBegin(GL_POLYGON);
        vec4 alphaColor = vec4(color, 0.0f);
        glColor4fv((GLfloat *) &alphaColor);
        glVertex3f(center.x - size/2, center.y - size/2, center.z+size/2);
        glVertex3f(center.x + size/2, center.y - size/2, center.z+size/2);
        glVertex3f(center.x + size/2, center.y + size/2, center.z+size/2);
        glVertex3f(center.x - size/2, center.y + size/2, center.z+size/2);
        glEnd();
    }
    
    void cube(vec3 center, float size){
        glPushMatrix();
        glTranslatef(center.x, center.y, center.z);
        //near
        vec4 alphaColor = vec4(blue, 0.0f);
        cubeSide(center, size, this->wallColors["near"]);
        glColor4fv((GLfloat*)&alphaColor);
        
        if (fractalsCompiled) {
            glCallList(KOCH_SNOWFLAKE);
        } else {
            glNewList(KOCH_SNOWFLAKE, GL_COMPILE);
                kochFlake();
            glEndList();
        }
        
        //right
        glRotatef(90, 0, 1, 0);
        cubeSide(center, size, this->wallColors["right"]);
        glColor4fv((GLfloat *)&red);
        
        if (fractalsCompiled) {
            glCallList(SQUARE_CARPET);
        } else {
            glNewList(SQUARE_CARPET, GL_COMPILE);
                squareCarpetFractal();
            glEndList();
        }
        
        //far
        glRotatef(90, 0, 1, 0);
        cubeSide(center, size, this->wallColors["far"]);
        
        if (fractalsCompiled) {
            glCallList(SIERPINSKY_TRIANGLE);
        } else {
            glNewList(SIERPINSKY_TRIANGLE, GL_COMPILE);
            sierpinskyTriangle();
            glEndList();
        }
        
        //left
        glRotatef(90, 0, 1, 0);
        cubeSide(center, size, this->wallColors["left"]);

        if (fractalsCompiled) {
            glCallList(HILBERT_CURVE);
        } else {
            glNewList(HILBERT_CURVE, GL_COMPILE);
            hilbertCurve();
            glEndList();
        }
        
        //bottom
        glRotatef(90, 1, 0, 0);
        cubeSide(center, size, this->wallColors["bottom"]);
        Fractal::chessFloor(center + vec3(-roomSize/2,-roomSize/2,roomSize/2 - roomSize/10000), roomSize, this->chessBoardTiles);
        
        
        
        //top
        glRotatef(180, 1, 0, 0);
        cubeSide(center, size, this->wallColors["top"]);
        sierpinskyCarpet();
        
        glPopMatrix();
        
        fractalsCompiled = true;
        
    }
    
    void kochFlake() {
        float flakePartSize = roomSize/6;
        glPushMatrix();
        glTranslatef(-roomSize/4, roomSize/4, roomSize/2 - roomSize/10000);
        Fractal::KochSnowFlake(4, flakePartSize);
        glTranslatef(flakePartSize*3, 0, 0);
        glRotatef(-90, 0, 0, 1);
        Fractal::KochSnowFlake(4, roomSize/6);
        glTranslatef(flakePartSize*3, 0, 0);
        glRotatef(-90, 0, 0, 1);
        Fractal::KochSnowFlake(4, roomSize/6);
        glTranslatef(flakePartSize*3, 0, 0);
        glRotatef(-90, 0, 0, 1);
        Fractal::KochSnowFlake(4, roomSize/6);
        glPopMatrix();

    }
    
    void squareCarpetFractal() {
        Fractal::SquareCarpet(center + vec3(-roomSize/2,-roomSize/2,roomSize/2-roomSize/10000), roomSize, 6);
    }
    
    void sierpinskyTriangle() {
        vec4 outerColor = vec4(azure, 0.0f);
        vec4 innerColor = vec4(black, 0.0f);
        glPushMatrix();
        
        glColor4fv((GLfloat *) &innerColor);
        
        for (int i = 0; i < 4; i++) {
            Fractal::SierpinskyLS(center + vec3(0,0,roomSize/2-roomSize/10000), roomSize/2, 5);
            glRotatef(90, 0, 0, 1);
        }
        
        glColor4fv((GLfloat *) &outerColor);
        for (int i = 0; i < 4; i++) {
            glPushMatrix();
            glRotatef(180, 0, 0, 1);
            glTranslatef(-roomSize/2, -roomSize/2, 0);
            Fractal::SierpinskyLS(center + vec3(0,0,roomSize/2-roomSize/10000), roomSize/2, 4);
            glPopMatrix();
            glRotatef(90, 0, 0, 1);
        }
        
        glPopMatrix();
    }
    
    void hilbertCurve() {
        //HILBERT CURVE
        glColor4fv((GLfloat *) &red);
        vec3 start = center + vec3(0,0,roomSize/2 - roomSize/10000);
        glLineWidth(roomSize/3);
        Fractal::HilbertCurve(start.x, start.y, start.z, start.x, roomSize, -roomSize, start.y, 5);
        glPushMatrix();
        glTranslatef(roomSize/2, -roomSize/2, 0);
        
        glBegin(GL_POINTS);
        glVertex3f(start.x, start.y, start.z);
        glEnd();
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i < hilbertPoints.size(); i++) {
            switch (i % 3) {
                case 0:
                    glColor4fv((GLfloat *) &red);
                    break;
                case 1:
                    glColor4fv((GLfloat *) &green);
                    break;
                case 2:
                    glColor4fv((GLfloat *) &magenta);
                    break;
                    
            }
            
            glVertex3f(hilbertPoints[i].x, hilbertPoints[i].y, hilbertPoints[i].z);
        }
        glEnd();
        glPopMatrix();
    }
    
    void sierpinskyCarpet() {
        vec4 alphaColor = vec4(yellow, 0.0f);
        glColor4fv((GLfloat *) &alphaColor);
        clock_t t = clock() - timer;
        double elapsed = (t/(double)CLOCKS_PER_SEC)*100;
        t = clock() - timer;
        elapsed = (t/(double)CLOCKS_PER_SEC)*100;
        Fractal::SierpinskySquareCarpet(center + vec3(-roomSize/2,-roomSize/2,roomSize/2 - roomSize/10000),roomSize,((int)elapsed%4)+1);
    }
    
};

class Camera {
public:
    Camera(vec3 eye = vec3(0,0,2), vec3 lookingAt = vec3(0,0,0), vec3 upSide = vec3(0,1,0)) {
        this->eye = eye;
        this->lookingAt = lookingAt;
        this->upSide = upSide;
        this->setDeltaMove(0.1);
        this->isAvatar = false;
        this->rotX = 0;
        this->rotY = 0;
        this->rotZ = 0;
        
        this->lookAt(eye, lookingAt, upSide);
        
    }
    void lookAt(vec3 eye, vec3 lookingAt, vec3 up) {
        gluLookAt(eye.x, eye.y, eye.z, lookingAt.x, lookingAt.y, lookingAt.z, up.x, up.y, up.z);
    }
    
    void setBall(Ball *ball) {
        this->ball = ball;
    }
    
    void setCoridor(Room *room) {
        this->room = room;
        vec3 center = room->getCenter();
    
        this->eye = center + vec3(0,0,room->getRoomSize()/2);
        this->lookingAt = center;
        this->upSide = vec3(0,1,0);
    }

    void move() {
        if (isAvatar) {
            vec3 eyeVert = ball->getPos() - vec3(10,10,10)*ball->getVel();
            vec3 ballVert = ball->getPos();
            vec3 normalizedMouseMove = vec3(normalizePosition(this->mouseX, this->mouseY),0.0f);
            
//            float roomSize = room->getRoomSize();
//            float angleX = normalizedMouseMove.y*-90;
//            float angleY = normalizedMouseMove.x*180;
//            float moveZ = abs(normalizedMouseMove.y)*roomSize/2;
//            printf("%f\n", moveZ);
////            rotate
//            mat4 rotateX = rotate(angleY, vec3(0.0f,1.0f,0.0f));
//            mat4 rotateY = rotate(angleX, vec3(1.0f,0.0f,0.0f));
//            mat4 finalRotMat =  rotateY*rotateX;
//            vec4 lookingAtVert = vec4(ballVert - vec3(0,0,roomSize/2),1.0f);
//            vec3 finalLookingAtVert = (vec3)(lookingAtVert*finalRotMat) +vec3(0,0,(-roomSize/2)+moveZ);
            
            
            this->lookAt(eyeVert, ballVert + normalizedMouseMove, upSide);

        } else {

            this->lookAt(eye, lookingAt, upSide);
            //rotate scene (camera moves only x-y) side jumps done by scene rotation
            glRotatef(90*this->rotX, 1, 0, 0);
            glRotatef(90*this->rotY, 0, 1, 0);
            glRotatef(90*this->rotZ, 0, 0, 1);
        }
        
    }
    
    void setAvatar(bool isAvatar) {
        this->isAvatar = isAvatar;
    }
    
    bool getCameraMode() {
        return this->isAvatar;
    }
    
    void setDeltaMove(float deltaMove) {
        this->deltaMove = deltaMove;
    }
    
    void setMouseMove(int x, int y) {
        this->mouseX = x;
        this->mouseY = y;
    }
    
    void moveUp() {
        
        float topBorder = room->getCenter().y + room->getRoomSize()/2;
        if ((eye.y + deltaMove) >= topBorder) {
            //negative modulo did not work properly
            if(rotY < 0) rotY = 4 - abs(rotY % 4);
            switch (rotY % 4) {
                case 0:
                    this->rotX++;
                    break;
                case 1:
                    this->rotZ++;
                    break;
                case 2:
                    this->rotX++;
                    break;
                case 3:
                    this->rotZ--;
                    break;
            }
            eye = eye * vec3(1,-1,1);
        } else {
            eye += vec3(0,deltaMove,0);
        }
        
    }
    
    void moveDown() {
        float bottomBorder = room->getCenter().y - room->getRoomSize()/2;
        if ((eye.y - deltaMove) <= bottomBorder) {
            if(rotY < 0) rotY = 4 - abs(rotY % 4);
            switch (rotY % 4) {
                case 0:
                    this->rotX--;
                    break;
                case 1:
                    this->rotZ--;
                    break;
                case 2:
                    this->rotX--;
                    break;
                case 3:
                    this->rotZ++;
                    break;
            }
            eye = eye * vec3(1,-1,1);
        } else {
            eye -= vec3(0,deltaMove,0);
        }
    }
    
    void moveLeft() {
        float leftBorder = room->getCenter().x - room->getRoomSize()/2;
        if ((eye.x - deltaMove) <= leftBorder) {
            if(rotX < 0) rotX = 4 - abs(rotX % 4);
            switch (rotX % 4) {
                case 0:
                    this->rotY++;
                    break;
                case 1:
                    this->rotZ--;
                    break;
                case 2:
                    this->rotY--;
                    break;
                case 3:
                    this->rotZ++;
                    break;
            }
            eye = eye * vec3(-1,1,1);
        } else {
            eye -= vec3(deltaMove,0,0);
        }
    }
    
    void moveRight() {
        float rightBorder = room->getCenter().x + room->getRoomSize()/2;
        if ((eye.x + deltaMove) >= rightBorder) {
            if(rotX < 0) rotX = 4 - abs(rotX % 4);
            switch (rotX % 4) {
                case 0:
                    this->rotY--;
                    break;
                case 1:
                    this->rotZ++;
                    break;
                case 2:
                    this->rotY++;
                    break;
                case 3:
                    this->rotZ--;
                    break;
            }
           
            eye = eye * vec3(-1,1,1);
        } else {
            eye += vec3(deltaMove,0,0);
        }
    }
    
    void moveForward() {
        eye -= vec3(0,0, deltaMove);
    }
    
    void moveBackward() {
        eye += vec3(0,0,deltaMove);
    }
private:
    vec3 eye, lookingAt, upSide;
    float deltaMove;
    bool isAvatar;
    Room *room;
    Ball *ball;
    int rotX, rotY, rotZ;
    int mouseX, mouseY;
};

class Scene {
public:
    Scene() {
        this->setProjection();
        
        camera = new Camera();
        
        room = new Room();
        room->setRoomSize(15.0);
        
        ball = new Ball();
        ball->setOuterBox(room->getRoomSize()/2);
        ball->setColor(blue);
        ball->setVel(vec3(0.04,0.05,-0.03));
        camera->setDeltaMove(0.4);
        camera->setCoridor(room);
        camera->setBall(ball);
    }
    
    
    
    void draw() {
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glEnable(GL_CULL_FACE);
        glCullFace(GL_FRONT);
        glBlendFunc(GL_ONE_MINUS_SRC_ALPHA,GL_SRC_ALPHA);
        glShadeModel(GL_FLAT);
        glClearColor(0.0, 0.0, 0.0, 0.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glLoadIdentity();
        camera->move();
        room->draw();
        ball->move();
        ball->draw();
        
        glutSwapBuffers();
    }
    
    void setProjection() {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glFrustum(-1, 1, -1, 1, 0.1, 20);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    }
    
    void keyPressed(unsigned char key) {
        switch (key) {
            //ESC
            case 27:
                exit(0);
                break;
            case 'w':
                camera->moveUp();
                break;
            case 's':
                camera->moveDown();
                break;
            case 'a':
                camera->moveLeft();
                break;
            case 'd':
                camera->moveRight();
                break;
            case 'f':
                camera->moveForward();
                break;
            case 'b':
                camera->moveBackward();
                break;
            case 'c':
                camera->setAvatar(!camera->getCameraMode());
                break;
            case '-':
                room->decrBoardcount();
                break;
            case '+':
                room->incrBoardCount();
                break;
            default:
                break;
        }
    }
    
    void mouseMove(int x, int y) {
        camera->setMouseMove(x, y);
    }
    
    void reshape(int w, int h) {
        glViewport (0, 0, (GLsizei) w, (GLsizei) h);
        glMatrixMode (GL_PROJECTION);
        glLoadIdentity ();
        gluPerspective(60.0, (GLfloat) w/(GLfloat) h, 0.1, 30.0);
        glMatrixMode (GL_MODELVIEW);
    }
    
private:
    Room *room;
    Camera *camera;
    Ball *ball;
};


Scene *scene;

void display(void) {
    scene->draw();
}

void keyboard(unsigned char key, int x, int y) {
    scene->keyPressed(key);
}

void menu() {
    printf("BALL LSD PRISON  Michal Gallovic\n\n");
    printf("Pohyb po stenach kocky: w/s/a/d\n");
    printf("Pohyb vpred a vzad: f/b (nesmiem byt v mode avatar)\n");
    printf("Prepinanie kamery: c\n");
    printf("Otacanie - pohyb mysou - ak som v mode avatar\n");
    printf("Pridavanie uberanie kachliciek na zemi: +/-\n");
    printf("\n\nENJOY");
}


void init() {
    scene = new Scene();
    menu();
}

void idle(void) {
    glutPostRedisplay();
}

void reshape(int w, int h) {
    scene->reshape(w, h);
}

void motion(int x, int y) {
    scene->mouseMove(x,y);
}


int main(int argc, char *argv[])
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowPosition(150, 150);
    glutInitWindowSize(600, 400);
    glutCreateWindow("Ball Prison");
    init();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);
    glutPassiveMotionFunc(motion);
    glutMainLoop();
    
    return 0;
}
