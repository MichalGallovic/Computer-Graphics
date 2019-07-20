#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <math.h>
#ifdef __APPLE__
#include "glm/glm.hpp"
#else
#include <glm/glm.hpp>
#endif


using namespace std;
using namespace glm;

struct Pixel {
    unsigned char R;
    unsigned char G;
    unsigned char B;
};

string imagePath;

class Camera {
public:
    Camera(vec3 eye = vec3(-1,2,2), vec3 lookingAt = vec3(0,0,0), vec3 upSide = vec3(0,1,0)) {
        this->eye = eye;
        this->lookingAt = lookingAt;
        this->upSide = upSide;
        this->setDeltaMove(0.1);
        this->lookAt(eye, lookingAt, upSide);
    }
    void lookAt(vec3 eye, vec3 lookingAt, vec3 up) {
        gluLookAt(eye.x, eye.y, eye.z, lookingAt.x, lookingAt.y, lookingAt.z, up.x, up.y, up.z);
    }

    void setDeltaMove(float deltaMove) {
        this->deltaMove = deltaMove;
    }
    
    void move() {
        this->lookAt(eye, lookingAt, upSide);
    }
    
    void moveUp() {
        eye += vec3(0,deltaMove,0);
    }
    
    void moveDown() {
        eye -= vec3(0,deltaMove,0);
    }
    
    void moveLeft() {
        eye -= vec3(deltaMove,0,0);
    }
    
    void moveRight() {
        eye += vec3(deltaMove,0,0);
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
};

class Texture{
public:
    Texture(int width, int height, vector<unsigned char> data,  GLuint texname)
    {
        this->width = width;
        this->height = height;
        texName = texname;
        glGenTextures(1, &texName);
        glBindTexture(GL_TEXTURE_2D, texName);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, 4, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, (unsigned char*)&data[0]);
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    }
    
    void renderTexture() {
        float aspect = 0.65;
        glBegin(GL_QUADS);
        glTexCoord2f(0.0f, 0.0f); glVertex3f(-1.0f, -1.0f*aspect,0);  // Bottom Left Of The Texture and Quad
        glTexCoord2f(1.0f, 0.0f); glVertex3f( 1.0f, -1.0f*aspect,0);  // Bottom Right Of The Texture and Quad
        glTexCoord2f(1.0f, 1.0f); glVertex3f( 1.0f,  1.0f*aspect,0);  // Top Right Of The Texture and Quad
        glTexCoord2f(0.0f, 1.0f); glVertex3f(-1.0f,  1.0f*aspect,0);  // Top Left Of The Texture and Quad
        glEnd();
    }
    
    void renderRoof() {
        glBegin(GL_QUADS);
        glTexCoord2f(0.0f, 0.0f); glVertex3f(-1.0f, -0.1f,0);  // Bottom Left Of The Texture and Quad
        glTexCoord2f(1.0f, 0.0f); glVertex3f( 1.0f, -0.1f,0);  // Bottom Right Of The Texture and Quad
        glTexCoord2f(1.0f, 1.0f); glVertex3f( 1.0f,  0.1f,0);  // Top Right Of The Texture and Quad
        glTexCoord2f(0.0f, 1.0f); glVertex3f(-1.0f,  0.1f,0);  // Top Left Of The Texture and Quad
        glEnd();
    }
    
    void renderPyramid() {
        glPushMatrix();
            glRotatef(180, 0, 0, 1);
            glRotatef(-45, 1, 0, 0);
            glTranslatef(0, 0, 1);
            renderTexture();
        glPopMatrix();
        
        glPushMatrix();
        glRotatef(90, 0, 1, 0);
        glRotatef(45, 0, 0, 1);
        glTranslatef(0, 0, 1);
        renderTexture();
        glPopMatrix();
    }

    int width,height;
    GLuint texName;
};

class ImageProcess {
public:
    ImageProcess(int width, int height, vector<unsigned char> rawImage) {
        this->width = width;
        this->height = height;
        this->image = this->RAW_RGB(rawImage);
    }
    void extractChannels(vector<unsigned char> tex) {
        
        for (int i = 0;  i< tex.size(); i+=3) {
            int index = i;
            R.push_back(tex[index]);
            R.push_back(0);
            R.push_back(0);
            
            G.push_back(0);
            G.push_back(tex[index+1]);
            G.push_back(0);
            
            B.push_back(0);
            B.push_back(0);
            B.push_back(tex[index+2]);
        }
        
    }
    
    vector<unsigned char> greyScaleRGB_RAW(vector<unsigned char> in){
        vector<unsigned char> out;
        for (int i = 0; i < in.size(); i+=3) {
            unsigned char Y = 0.2126*in[i] + 0.7152*in[i+1] + 0.0722*in[i+2];
            out.push_back(Y);
            out.push_back(Y);
            out.push_back(Y);
        }
        
        return out;
    }
    
    vector<vector<Pixel>> RAW_RGB(vector<unsigned char> in) {
        vector<vector<Pixel>> image;
        for (int row = 0; row < height; row++) {
            vector<Pixel> image_row;
            for (int col = 0; col < width; col++) {
                int index = row*width*3 + col*3;
                struct Pixel pix;
                pix.R = in[index];
                pix.G = in[index+1];
                pix.B = in[index+2];
                image_row.push_back(pix);
            }
            image.push_back(image_row);
        }
        return image;
    }
    
    vector<vector<Pixel>> tresholdSky(vector<vector<Pixel>> in) {
        vector<vector<Pixel>> image;
        vec3 average;
        average.x = 0;
        average.y = 0;
        average.z = 0;
        int row = 0, col = 0;
        int counter = 0;

        // histogram
        
        if (height > 200) {
            //pyramid
            for (row = height-100; row >  height-200; row--) {
                for (col = 100; col < 200; col++) {
                    average.x += in[row][col].R;
                    average.y += in[row][col].G;
                    average.z += in[row][col].B;
                    counter++;
                }
            }
        } else {
            // rooftop
            for (row = height-1; row >  height-60; row--) {
                for (col = 300; col < 400; col++) {
                    average.x += in[row][col].R;
                    average.y += in[row][col].G;
                    average.z += in[row][col].B;
                    counter++;
                }
            }
        }
        
        Pixel averagePix;
        averagePix.B = average.z/counter;
        
        
        for(int row = 0; row < height; row++) {
            vector<Pixel> image_row;
            for (int col = 0; col < width; col++) {
                int index = row*width*3 + col*3;
                Pixel pix = in[row][col];
                if (pix.B > averagePix.B || pix.G < 90) {
                    pix.R = 0;
                    pix.G = 0;
                    pix.B = 0;
                }
                image_row.push_back(pix);
                
            }
            image.push_back(image_row);
        }
        
        return image;
        
    }
    
    vector<vector<Pixel>> eraseAreas(vector<vector<Pixel>> in) {
        vector<vector<Pixel>> image;
        
        for(int row = 0; row < height; row++) {
            vector<Pixel> image_row;
            for (int col = 0; col < width; col++) {

                Pixel pix = in[row][col];
                if (col < 120 && row > 80 && row < 100) {
                    pix.R = 0;
                    pix.G = 0;
                    pix.B = 0;
                }
                
                if (col < 75 && row > 50 && row < 100) {
                    pix.R = 0;
                    pix.G = 0;
                    pix.B = 0;
                }
                
                if (col > 830 && row > 70 && row < 100) {
                    pix.R = 0;
                    pix.G = 0;
                    pix.B = 0;
                }
                if (col > 780 && row > 100 && row < 120) {
                    pix.R = 0;
                    pix.G = 0;
                    pix.B = 0;
                }
                
                image_row.push_back(pix);
                
            }
            image.push_back(image_row);
        }
        
        return image;
    }
    
    vector<unsigned char> RGB_RAW(vector<vector<Pixel>> in) {
        vector<unsigned char> out;
        for (int row = 0; row < height; row++) {
            for (int col = 0; col < width; col++) {
                Pixel pix = in[row][col];
                
                out.push_back(pix.R);
                out.push_back(pix.G);
                out.push_back(pix.B);
                if ((pix.R == 0 && pix.G == 0 && pix.B == 0)) {
                    out.push_back(0);
                } else {
                    if (row >= 1) {
                        out.push_back(255);
                    } else {
                        out.push_back(0);
                    }
                }
                

                
            }
        }
        

        return out;
    }
    
    vector<vector<Pixel>> cropImage(vector<vector<Pixel>> in, int leftBorder, int rightBorder, int bottomBorder, int topBorder) {
        vector<vector<Pixel>> image;
        for(int row = 0; row < height; row++) {
            if (row > bottomBorder && row < topBorder) {
                vector<Pixel> img_row;
                vector<Pixel> original_image_row;
                for (int col = 0; col < width; col++) {
                    Pixel pix = in[row][col];
                    if (col > leftBorder && col < rightBorder) {
                        img_row.push_back(pix);
                        original_image_row.push_back(this->image[row][col]);
                    }
                }
                image.push_back(img_row);
                this->croppedImage.push_back(original_image_row);
            }
        }
        
        this->width = rightBorder - leftBorder - 1;
        this->height = topBorder - bottomBorder - 1;
        
        return image;
    }
    
    vector<vector<Pixel>> fillEmptySpaces(vector<vector<Pixel>> in){
        // find center position of kernel (half of kernel size)
        int kCols = 7;
        int kRows = 7;
        int kCenterX = kCols / 2;
        int kCenterY = kRows / 2;
        int rows = height;
        int cols = width;
        int mm,nn, ii,jj;
        float kernel[7][7] = {
            {1,1,1,1,1,1,1},
            {1,1,1,1,1,1,1},
            {1,1,1,1,1,1,1},
            {1,1,1,1,1,1,1},
            {1,1,1,1,1,1,1},
            {1,1,1,1,1,1,1},
            {1,1,1,1,1,1,1}
        };
        vector<vector<Pixel>> image = in;
        vector<vector<vec3>> out (height, vector<vec3>(width));
        
        for(int i=0; i < rows; ++i)              // rows
        {
            for(int j=0; j < cols; ++j)          // columns
            {
                
                for(int m=0; m < kRows; ++m)     // kernel rows
                {
                    mm = kRows - 1 - m;      // row index of flipped kernel
                    
                    for(int n=0; n < kCols; ++n) // kernel columns
                    {
                        nn = kCols - 1 - n;  // column index of flipped kernel
                        
                        // index of input signal, used for checking boundary
                        ii = i + (m - kCenterY);
                        jj = j + (n - kCenterX);
                        
                        // ignore input samples which are out of bound
                        if( ii >= 0 && ii < rows && jj >= 0 && jj < cols ) {
                            out[i][j].x += image[ii][jj].R * kernel[mm][nn];
                            out[i][j].y += image[ii][jj].G * kernel[mm][nn];
                            out[i][j].z += image[ii][jj].B * kernel[mm][nn];
                            
                        }
                    }
                }
                if (out[i][j].x == image[i][j].R) {
                    out[i][j].x = 0;
                } else {
                    out[i][j].x = this->croppedImage[i][j].R;
                }
                if (out[i][j].y == image[i][j].G) {
                    out[i][j].y = 0;
                } else {
                    out[i][j].y = this->croppedImage[i][j].G;
                }
                if (out[i][j].z == image[i][j].B) {
                    out[i][j].z = 0;
                } else {
                    out[i][j].z = this->croppedImage[i][j].B;
                }
            }
        }
        
        
        
        return this->vecToRGB(out);
    }
    
    vector<vector<Pixel>> outlinePyramid(vector<vector<Pixel>> in) {
        vector<vector<Pixel>> image;
        
        for(int row = 0; row < height; row++) {
            vector<Pixel> image_row;
            for (int col = 0; col < width; col++) {
                Pixel pix = in[row][col];
                if ((pix.B > 180)) {
                    
                    pix.R = 255;
                    pix.G = 0;
                    pix.B = 0;
                }
                image_row.push_back(pix);
                
            }
            image.push_back(image_row);
        }
        
        return image;
    }
    
    vector<vector<Pixel>> outlinePyramid2(vector<vector<Pixel>> in) {
        vector<vector<Pixel>> image;
        
        for(int row = 0; row < height; row++) {
            vector<Pixel> image_row;
            for (int col = 0; col < width; col++) {
                Pixel pix = in[row][col];
                int treshold = 30;
                if ((pix.B > treshold && pix.G > treshold && pix.B > treshold)) {
                    
                    pix.R = 255;
                    pix.G = 0;
                    pix.B = 0;
                }
                image_row.push_back(pix);
                
            }
            image.push_back(image_row);
        }
        
        return image;
    }
    
    vector<vector<Pixel>> process(int size, float kernel[3][3], vector<vector<Pixel>> in){
        // find center position of kernel (half of kernel size)
        int kCols = size;
        int kRows = size;
        int kCenterX = kCols / 2;
        int kCenterY = kRows / 2;
        int rows = height;
        int cols = width;
        int mm,nn, ii,jj;
        vector<vector<Pixel>> image = in;
        vector<vector<vec3>> out (height, vector<vec3>(width));
        Pixel pix;
        for(int i=0; i < rows; ++i)              // rows
        {
            for(int j=0; j < cols; ++j)          // columns
            {
                
                for(int m=0; m < kRows; ++m)     // kernel rows
                {
                    mm = kRows - 1 - m;      // row index of flipped kernel
                    
                    for(int n=0; n < kCols; ++n) // kernel columns
                    {
                        nn = kCols - 1 - n;  // column index of flipped kernel
                        
                        // index of input signal, used for checking boundary
                        ii = i + (m - kCenterY);
                        jj = j + (n - kCenterX);
                        
                        // ignore input samples which are out of bound
                        if( ii >= 0 && ii < rows && jj >= 0 && jj < cols ) {
                            out[i][j].x += image[ii][jj].R * kernel[mm][nn];
                            out[i][j].y += image[ii][jj].G * kernel[mm][nn];
                            out[i][j].z += image[ii][jj].B * kernel[mm][nn];
                            
                        }
                    }
                }
                
            }
        }
        
        apply_convolution_rules(&out);
        
        return vecToRGB(out);
    }
    
    float angleFromSobel(vector<vector<Pixel>> in) {
        float sobel_hor [3][3]= {
            {1,2,1},
            {0,0,0},
            {-1,-2,-1}
        };
        int height = (int)in.size();
        int width = (int)in[0].size();
        
        float sobel_vert [3][3] = {
            {1,0,-1},
            {2,0,-2},
            {1,0,-1}
        };
        float sum_hor = 0, sum_vert = 0 ,val = 0;
        int kCols = 3;
        int kRows = 3;
        int kCenterX = kCols / 2;
        int kCenterY = kRows / 2;
        int rows = height;
        int cols = width;
        int mm,nn, ii,jj;
        vector<vector<Pixel>> image = in;
        vector<vector<vec3>> out_hor (height, vector<vec3>(width));
        vector<vector<vec3>> out_vert (height, vector<vec3>(width));
        
        for(int i=0; i < rows; ++i)              // rows
        {
            for(int j=0; j < cols; ++j)          // columns
            {
                
                for(int m=0; m < kRows; ++m)     // kernel rows
                {
                    mm = kRows - 1 - m;      // row index of flipped kernel
                    
                    for(int n=0; n < kCols; ++n) // kernel columns
                    {
                        nn = kCols - 1 - n;  // column index of flipped kernel
                        
                        // index of input signal, used for checking boundary
                        ii = i + (m - kCenterY);
                        jj = j + (n - kCenterX);
                        
                        // ignore input samples which are out of bound
                        if( ii >= 0 && ii < rows && jj >= 0 && jj < cols ) {
                            out_hor[i][j].x += image[ii][jj].R * sobel_hor[mm][nn];
                            out_hor[i][j].y += image[ii][jj].G * sobel_hor[mm][nn];
                            out_hor[i][j].z += image[ii][jj].B * sobel_hor[mm][nn];
                            
                            out_vert[i][j].x += image[ii][jj].R * sobel_vert[mm][nn];
                            out_vert[i][j].y += image[ii][jj].G * sobel_vert[mm][nn];
                            out_vert[i][j].z += image[ii][jj].B * sobel_vert[mm][nn];
                            
                        }
                    }
                }
                val = ceil(std::sqrt((out_hor[i][j].x * out_hor[i][j].x) + (out_vert[i][j].x * out_vert[i][j].x)));
                if (val > 50) {
                    sum_hor+= out_hor[i][j].x;
                    sum_vert+= out_vert[i][j].x;
                }
                
            }
        }
        

        return abs((atan(sum_vert/sum_hor) *180) / M_PI);
        
    }
    
    void apply_convolution_rules(vector<vector<vec3>> *pix_after_convol) {
        for (int row = 0; row < height; row++) {
            for (int col = 0; col < width; col++) {
                if ((*pix_after_convol)[row][col].x > 255) {
                    (*pix_after_convol)[row][col].x = 255;
                }
                if ((*pix_after_convol)[row][col].y > 255) {
                    (*pix_after_convol)[row][col].y = 255;
                }
                if ((*pix_after_convol)[row][col].z > 255) {
                    (*pix_after_convol)[row][col].z = 255;
                }
                
                if ((*pix_after_convol)[row][col].x < 0) {
                    (*pix_after_convol)[row][col].x = 0;
                }
                if ((*pix_after_convol)[row][col].y < 0) {
                    (*pix_after_convol)[row][col].y = 0;
                }
                if ((*pix_after_convol)[row][col].z < 0) {
                    (*pix_after_convol)[row][col].z = 0;
                }
            }
        }
        
    }
    
    int width, height;
    vector<unsigned char> R,G,B;
    vector<vector<Pixel>> image, croppedImage, processed;
private:
    vector<vector<Pixel>> vecToRGB(vector<vector<vec3>> in) {
        vector<vector<Pixel>> out (height, vector<Pixel>(width));
        for(int row = 0; row < height; row++) {
            for (int col = 0; col < width; col++) {
                Pixel pix;
                pix.R = in[row][col].x;
                pix.G = in[row][col].y;
                pix.B = in[row][col].z;
                out[row][col] = pix;
            }
        }
        return out;
    }
};

class Pyramid {
public:
    Pyramid(float angle, int textureAngle, float distFromCenter, Texture *tex) {
        this->angle = angle;
        this->textureAngle = textureAngle;
        this->distFromCenter = distFromCenter;
        this->texture = tex;
    }
    
    void draw() {
        glBindTexture(GL_TEXTURE_2D, texture->texName);
        // pyramida
        glPushMatrix();
        for (int i = 0; i < 4; i++) {
            glPushMatrix();
            glTranslatef(0, 0, distFromCenter);
            glRotatef(-angle, 1, 0, 0);
            glRotatef(-textureAngle, 0, 0, 1);
            texture->renderTexture();
            glPopMatrix();
            glRotatef(90, 0, 1, 0);
        }
        glPopMatrix();
        
    }
    float angle;
    float distFromCenter;
    int textureAngle;
    Texture *texture;
};

class BMP{
public:
    vector<unsigned char> readBMP(char* filename)
    {
        int i;
        FILE* f = fopen(filename, "rb");
        unsigned char info[54];
        fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header
        
        // extract image height and width from header
        this->width = abs(*(int*)&info[18]+1);
        this->height = abs(*(int*)&info[22]);
        
        int size = 3 * width * height;
        unsigned char* data = new unsigned char[size]; // allocate 3 bytes per pixel
        fread(data, sizeof(unsigned char), size, f); // read the rest of the data at once
        fclose(f);
        
        for(i = 0; i < size; i += 3)
        {
            unsigned char tmp = data[i];
            data[i] = data[i+2];
            data[i+2] = tmp;
            rawRGB.push_back(data[i]);
            rawRGB.push_back(data[i+1]);
            rawRGB.push_back(data[i+2]);
        }
        
        return rawRGB;
    }
    
    
    vector<unsigned char> rawRGB;
    int width,height;
};

class Scene {
public:
    Scene() {
        this->setProjection();
        this->settingsEnable();
        this->counter = 0.02;
        this->counter2 = 0.68;
        float sobel_hor [3][3]= {
            {1,2,1},
            {0,0,0},
            {-1,-2,-1}
        };
        
        float sobel_vert [3][3] = {
            {1,0,-1},
            {2,0,-2},
            {1,0,-1}
        };
        
        float laplace [3][3] = {
            {0,-1,0},
            {-1,4,-1},
            {0,-1,0}
        };

        char* imageFilePath = &imagePath[0u];
        
        // load image + parse image
        imageFile = new BMP();
        imageFile->readBMP(imageFilePath);
        
        // image processing
        image = new ImageProcess(imageFile->width, imageFile->height, imageFile->rawRGB);
        image->processed = image->fillEmptySpaces(image->eraseAreas(image->tresholdSky(image->cropImage(image->image, 115, 1012, 50, 300))));
//        image->processed = image->fillEmptySpaces(image->tresholdSky(image->cropImage(image->image, 115, 1012, 50, 300)));
        image->processed = image->outlinePyramid(image->processed);

        
        // texture init
        texture = new Texture(image->width, image->height, image->RGB_RAW(image->processed),1);

        // angle computing
        BMP *imageFile2 = new BMP();
        imageFile2->readBMP(imageFilePath);
        ImageProcess *image2 = new ImageProcess(imageFile2->width, imageFile2->height, imageFile2->rawRGB);
        image2->processed = image2->RAW_RGB(image2->greyScaleRGB_RAW(imageFile2->rawRGB));
        image2->processed = image2->cropImage(image2->processed, 785, 845, 408, 468);
        float angleFix = 20.0f;
        float computedAngle = image2->angleFromSobel(image2->processed) + angleFix;
        
        // pyramid
        pyramid = new Pyramid(computedAngle, 2, 0.65, texture);
        
        ImageProcess *image3 = new ImageProcess(imageFile->width, imageFile->height, imageFile->rawRGB);
        image3->processed = image3->fillEmptySpaces(image3->eraseAreas(image3->tresholdSky(image3->cropImage(image3->image, 115, 1012, 290, 350))));
        image3->processed = image3->outlinePyramid(image3->processed);
        texture2 = new Texture(image3->width, image3->height, image3->RGB_RAW(image3->processed),2);
        // camera
        camera = new Camera();
//
        
        

    }
    
    void settingsEnable() {
        glEnable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_TEXTURE_2D);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable( GL_ALPHA_TEST );
        glAlphaFunc( GL_GREATER, 0 );
    }
    
    void settingsDisable() {
        glDisable(GL_BLEND);
        glDisable(GL_ALPHA_TEST);
        glDisable(GL_TEXTURE_2D);
    }
    
    void draw() {
        glClearColor(0.0, 0.0, 0.0, 1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glLoadIdentity();

        camera->move();
        
        pyramid->draw();
        
        glBindTexture(GL_TEXTURE_2D, 2);
        
        glPushMatrix();
        glTranslatef(0, 0.54, 0);
        for (int i = 0; i < 4; i++) {

            glPushMatrix();
            glTranslatef(0.024, 0, 0.2);
            glScalef(0.98, 0.98, 0.98);
            glRotatef(-2.1, 0, 0, 1);
            texture2->renderRoof();
            glPopMatrix();
            glRotatef(90, 0, 1, 0);
            
        }
        glPopMatrix();
        
        this->settingsDisable();
        GLfloat size = 0.2;
        
        glPushMatrix();
        glTranslatef(0, 0.63, 0);
        
        glBegin(GL_POLYGON);
            glColor3f(0.4, 0.4, 0.4);
        glVertex3f(-size, 0, size);
        glVertex3f(size, 0, size);
        glVertex3f(size, 0, -size);
        glVertex3f(-size, 0, -size);
        glEnd();
        glPopMatrix();

        glPushMatrix();
        glColor3f(1.0, 0.4, 0.4);
        
        glTranslatef(0.0, 0.54, 0.0);
        for (int i = 0; i < 4; i++) {
            glPushMatrix();
            glTranslatef(0.15, 0, 0.145);
                glScalef(0.7, 1.7, 1.0);
                glutSolidCube(0.1);
            glPopMatrix();
            glRotatef(90, 0, 1, 0);
        }
        glPopMatrix();
        
        glPushMatrix();
        glTranslatef(0, -0.35, 0);
        glColor3f(0.3, 0.3, 0.3);
        glBegin(GL_POLYGON);
        glVertex3f(-10, 0, -10);
        glVertex3f(10, 0, -10);
        glVertex3f(10, 0, 10);
        glVertex3f(-10, 0, 10);
        glEnd();
        glPopMatrix();
        
        
        this->settingsEnable();
        
       
        


        

        
        
        
        
        glutSwapBuffers();
    }
    
    void reshape(int w, int h) {
        glViewport (0, 0, (GLsizei) w, (GLsizei) h);
        glMatrixMode (GL_PROJECTION);
        glLoadIdentity ();
        
        gluPerspective(60.0, (GLfloat) w/(GLfloat) h, 0.1, 30.0);
        glMatrixMode (GL_MODELVIEW);
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
            case '+':
                this->counter+=0.001;
                 printf("%f ", this->counter);
                break;
            case '-':
                this->counter-=0.001;
                 printf("%f ", this->counter);
                break;
            case 'k':
                this->counter2+=0.01;
                 printf("%f ", this->counter2);
                break;
            case 'l':
                this->counter2-=0.01;
                 printf("%f ", this->counter2);
                break;
            case 'f':
                camera->moveForward();
                break;
            case 'b':
                camera->moveBackward();
                break;
        }
       
    }
    
    void setProjection() {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glFrustum(-1, 1, -1, 1, 0.1, 20);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    }
    
    BMP *imageFile;
    ImageProcess *image;
    Texture *texture;
    Texture *texture2;
    Camera *camera;
    Pyramid *pyramid;
    float counter;
    float counter2;

};

Scene *scene;

void display(void) {
    scene->draw();
}

void init() {
    scene = new Scene();
}

void repeat() {
    glutPostRedisplay();
}

void reshape(int w, int h) {
    scene->reshape(w, h);
}

void keyboard(unsigned char key, int x, int y) {
    scene->keyPressed(key);
}

int main(int argc, char **argv)
{
    if (argc > 1) {
        imagePath = argv[1];
    } else {
        printf("Please provide absolute path to image");
        return 0;
    }
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutInitWindowSize(1000, 600);
    glutCreateWindow("Image processing");
    init();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(repeat);
    glutMainLoop();

    return 0;
}





