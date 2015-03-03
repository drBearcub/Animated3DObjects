//David Chen
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <OpenGL/gl.h>
 #include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

#define PI 3.14159265

using namespace Eigen;
using namespace std;

void lineDDA (int x0, int y0, int xEnd, int yEnd, float* board, float r0, float g0, float b0, float rEnd, float gEnd, float bEnd);
void makePix(int x, int y, float* board, float r, float g, float b) ;
void rasterize(int y, float* board);

float Ka,Kd,Ks,IA,IL,K,n;
int windowSize;
float* mainBuffer;
int objectNum = 0;
int lightSource[3];
int viewingVector_XY[3];
int viewingVector_YZ[3];
int viewingVector_XZ[3];

class Vertex{
public:
  float x,y,z;
  float r,g,b;
  float normal[3];

  Vertex(){
    x = 0;
    y = 0;
    z = 0;
  }

  Vertex(float x, float y, float z){
    this->x = x;
    this->y = y;
    this->z = z;
  }
};

class SurfacePoints{
public:
  int p1,p2,p3;
  SurfacePoints(int p1, int p2, int p3){
    this->p1 = p1;
    this->p2 = p2;
    this->p3 = p3;
  }
};

class Surface {
public:
  vector<Vertex*> points;

  float minX = 900; float minY = 900; float minZ = 900;
  float maxX = 0; float maxY = 0; float maxZ = 0;
  float avgX = 0; float avgY = 0; float avgZ = 0;
  float normal[3];

  int p1,p2,p3;

  Surface(Vertex* a, Vertex* b, Vertex* c){
    points.push_back(a);
    points.push_back(b);
    points.push_back(c);    
    computeValues();
  }


  void computeValues(){
    for(int i = 0; i < 3; ++i){
      minX = min(minX, points[i]->x);
      minY = min(minY, points[i]->y);
      minZ = min(minZ, points[i]->z);  
      maxX = max(maxX, points[i]->x);
      maxY = max(maxY, points[i]->y);
      maxZ = max(maxZ, points[i]->z);

      float a[3];
      float b[3];

      a[0] = (points[0]->x - points[1]->x);
      a[1] = (points[0]->y - points[1]->y);
      a[2] = (points[0]->z - points[1]->z);
      b[0] = (points[2]->x - points[1]->x);
      b[1] = (points[2]->y - points[1]->y);
      b[2] = (points[2]->z - points[1]->z);

      normal[0] = a[1]*b[2] - a[2]*b[1];
      normal[1] = a[2]*b[0] - a[0]*b[2];
      normal[2] = a[0]*b[1] - a[1]*b[0];  

      avgX += points[i]->x;
      avgY += points[i]->y;
      avgZ += points[i]->z;      
    }
    
  }
};

class Object {
public:
  vector<Surface> surfaces;
  vector<Vertex> points;
  float normals[4][3];
  vector<SurfacePoints> surfpoints;

  float xCenter = 0; float yCenter = 0; float zCenter = 0;
  float xMin = 10;
  float yMin = 10;
  float zMin = 10;

  Object(){
  }

  void computeValues(){
    for(int i = 0; i < points.size(); ++i) {
      xMin = min(points[i].x,xMin);
      yMin = min(points[i].x,yMin);
      zMin = min(points[i].x,zMin);  
      xCenter += points[i].x;
      yCenter += points[i].y;
      zCenter += points[i].z;      
    }
    xCenter = xCenter/4;
    yCenter = yCenter/4;
    zCenter = zCenter/4;
  }

  void addVertex(float x, float y, float z) {
      Vertex temp(x,y,z);
      points.push_back(temp);
  }

  void addSurface(int one, int two, int three){
    SurfacePoints temp(one,two,three);
    surfpoints.push_back(temp);
  }

  void generateSurfaces() {
    for(int i = 0; i < surfpoints.size(); ++i) {
      Surface temp(
        &points[surfpoints[i].p1-1],
        &points[surfpoints[i].p2-1],
        &points[surfpoints[i].p3-1]);

      temp.p1 = surfpoints[i].p1;
      temp.p2 = surfpoints[i].p2;
      temp.p3 = surfpoints[i].p3;

      surfaces.push_back(temp);
    }
  }

  void computeVertexNormal(){ //loop through the surfaces and add values to normal, then average out by 3.
    for(int i = 0; i < surfaces.size(); ++i) {
      points[surfaces[i].p1-1].normal[0] += surfaces[i].normal[0]; //x
      points[surfaces[i].p1-1].normal[1] += surfaces[i].normal[1]; //y
      points[surfaces[i].p1-1].normal[2] += surfaces[i].normal[2]; //z 
           
      points[surfaces[i].p2-1].normal[0] += surfaces[i].normal[0];
      points[surfaces[i].p2-1].normal[1] += surfaces[i].normal[1];
      points[surfaces[i].p2-1].normal[2] += surfaces[i].normal[2]; 

      points[surfaces[i].p3-1].normal[0] += surfaces[i].normal[0];
      points[surfaces[i].p3-1].normal[1] += surfaces[i].normal[1];
      points[surfaces[i].p3-1].normal[2] += surfaces[i].normal[2]; 
    }


    for(int i = 0; i < 4; ++i) {
      float length;
      points[i].normal[0] = points[i].normal[0]/3;
      points[i].normal[1] = points[i].normal[1]/3;
      points[i].normal[2] = points[i].normal[2]/3; 

      length = sqrt(points[i].normal[0]*points[i].normal[0] +
                    points[i].normal[1]*points[i].normal[1] +
                    points[i].normal[2]*points[i].normal[2]);

      points[i].normal[0] = points[i].normal[0]/length;
      points[i].normal[1] = points[i].normal[1]/length;
      points[i].normal[2] = points[i].normal[2]/length;
    }

  }

  void computeVertexIntensity(){ // TO DO
    float l[3];
    float diffusedMultiplier;
    float length;
    for(int i = 0; i < points.size(); ++i) {
      l[0] = lightSource[0]-points[i].x;
      l[1] = lightSource[1]-points[i].y;
      l[2] = lightSource[2]-points[i].z;   
      length = sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
      l[0] = l[0]/length;
      l[1] = l[1]/length;
      l[2] = l[2]/length; 
      diffusedMultiplier = l[0] * points[i].normal[0] + l[1]* points[i].normal[1] + l[2] * points[i].normal[2];
      diffusedMultiplier *= IL*Kd * 2;
      points[i].r = 0.2 + diffusedMultiplier;
      points[i].g = 0.2 + diffusedMultiplier;
      points[i].b = 0.4 + diffusedMultiplier;
    }

  }

  void translation(float dx, float dy, float dz) {
    for(int i = 0; i < points.size(); ++i) {
      points[i].x += dx;
      points[i].y += dy;
      points[i].z += dz;        
    }
  }

  void scale(float s) {
    for(int i = 0; i < points.size(); ++i) {
      points[i].x = (points[i].x - xCenter) * s + xCenter;
      points[i].y = (points[i].y - yCenter) * s + yCenter;
      points[i].z = (points[i].z - zCenter) * s + zCenter;        
    }
  }

    void rotationUpInto(float alpha, int direction){ // 1-up, 2-down,3-left,4-right
      alpha = (alpha * PI)/180;
        if(direction == 1){
          translation(-xCenter,-yCenter,-zCenter);
          MatrixXf Rx(4,4);
          Rx << 1 , 0 , 0 , 0
               , 0 , cos(alpha), -1*sin(alpha) , 0
               , 0 , sin(alpha), cos(alpha)    , 0
               , 0 , 0 , 0 , 1;
          MatrixXf result(4,1);
          for(int i = 0; i < points.size(); ++i) {
            MatrixXf Pt(4,1);
            Pt << points[i].x , points[i].y, points[i].z , 1;
            result = Rx * Pt; 
            points[i].x = result(0,0);
            points[i].y = result(1,0);
            points[i].z = result(2,0);
          } 
        }else if(direction == 3){
          translation(-xCenter,-yCenter,-zCenter);
          MatrixXf Rx(4,4);
          Rx <<  cos(alpha) , 0 , -1*sin(alpha) , 0
               , 0 , 1, 0, 0
               , sin(alpha), 0, cos(alpha), 0
               , 0 , 0 , 0 , 1;
          MatrixXf result(4,1);
          for(int i = 0; i < points.size(); ++i) {
            MatrixXf Pt(4,1);
            Pt << points[i].x , points[i].y, points[i].z , 1;
            result = Rx * Pt;
            points[i].x = result(0,0);
            points[i].y = result(1,0);
            points[i].z = result(2,0);
          } 
        }
        translation(xCenter, yCenter, zCenter);
    }

};

vector<Object> Objects;
vector<Surface*> Triangles;

class Scene {
public:
  float normalizeRatio;
    float normalizeSize;
    float transferred = false;

  Scene(float normalizeSize = 0.5) {
    this->normalizeSize = normalizeSize;
  }

  void allTrianglesRecompute() {
    for(int i = 0; i < Triangles.size(); ++i) {
      Triangles[i]->computeValues();
    }
  }

  void allObjectsRecompute() {
    for(int i = 0; i < Objects.size(); ++i) {
      Objects[i].computeValues();
    }
  }

  void transferTriangles(){
    if(transferred == false){ // protection against duplication in Triangles array
      for(int i = 0; i < Objects.size(); ++i) {
        for(int j = 0; j < Objects[i].surfaces.size(); ++j) {
          Triangles.push_back(&Objects[i].surfaces[j]);
        }
      }
      transferred = true;
    }
  }


  void allVertexRecompute(){
    for(int i = 0; i < Objects.size(); ++i) {
      Objects[i].computeVertexNormal();
      Objects[i].computeVertexIntensity();
    }
  }


  void configurationDisplayed() {
    transferred = false;
    Triangles.clear();
  }

  void sortSurfaces(int order){

    vector<Surface*> temp;
    if(order == 1) { //XY, sort Z
      for(int j = 0; j < Triangles.size(); ++j) {
        if(temp.size() == 0){
          temp.push_back(Triangles[j]);
        } else if(temp[temp.size()-1]->minZ <= Triangles[j]->minZ) {
          temp.push_back(Triangles[j]);
        } else if(temp[0]->minZ >= Triangles[j]->minZ) {
          temp.insert(temp.begin(),Triangles[j]);
        } else {
          for(int k = 0; k < temp.size()-1; k++) {
            if(temp[k]->minZ <= Triangles[j]->minZ && Triangles[j]->minZ <= temp[k+1]->minZ ) {
              temp.insert(temp.begin() + k, Triangles[j]);
              break;
            }
          }
        }
      }
      Triangles = temp;
    }

    if(order == 2) { 
      for(int j = 0; j < Triangles.size(); ++j) {
        if(temp.size() == 0){
          temp.push_back(Triangles[j]);
        } else if(temp[temp.size()-1]->minY <= Triangles[j]->minY) {
          temp.push_back(Triangles[j]);
        } else if(temp[0]->minY >= Triangles[j]->minZ) {
          temp.insert(temp.begin(),Triangles[j]);
        } else {
          for(int k = 0; k < temp.size()-1; k++) {
            if(temp[k]->minY <= Triangles[j]->minY && Triangles[j]->minY <= temp[k+1]->minY ) {
              temp.insert(temp.begin() + k, Triangles[j]);
              break;
            }
          }
        }
      }
      Triangles = temp;
    }

    if(order == 3) { //XY, sort Z
      for(int j = 0; j < Triangles.size(); ++j) {
        if(temp.size() == 0){
          temp.push_back(Triangles[j]);
        } else if(temp[temp.size()-1]->minX <= Triangles[j]->minX) {
          temp.push_back(Triangles[j]);
        } else if(temp[0]->minX >= Triangles[j]->minX) {
          temp.insert(temp.begin(),Triangles[j]);
        } else {
          for(int k = 0; k < temp.size()-1; k++) {
            if(temp[k]->minX <= Triangles[j]->minX && Triangles[j]->minX <= temp[k+1]->minX ) {
              temp.insert(temp.begin() + k, Triangles[j]);
              break;
            }
          }
        }
      }
      Triangles = temp;
    }

  }

  void normalize_XY() {
    float xMin = 9999; float xMax = -9999; float yMin = 9999; float yMax = -9999;
    for(int i = 0; i < Triangles.size(); i++) {
      for(int j = 0; j < Triangles[i]->points.size(); ++j) {
        xMin = min(xMin, Triangles[i]->points[j]->x);
        xMax = max(xMax, Triangles[i]->points[j]->x);
        yMin = min(yMin, Triangles[i]->points[j]->y);
        yMax = max(yMax, Triangles[i]->points[j]->y);               
      }   
    }
    float deltaX = xMax - xMin;
    float deltaY = yMax - yMin;
    normalizeRatio = normalizeSize/max(deltaY,deltaX);

    for(int i = 0; i < Triangles.size(); ++i) {
      if(Triangles[i]->normal[2] <=  0) {
        float* tempBuffer = new float [windowSize * windowSize * 3 * 4];

        for(int j = 0; j < 2; ++j){
          lineDDA( 
                  (Triangles[i]->points[j]->x - xMin) * normalizeRatio * windowSize, 
                  (Triangles[i]->points[j]->y - yMin) * normalizeRatio * windowSize, 
                  (Triangles[i]->points[j+1]->x - xMin) * normalizeRatio * windowSize, 
                  (Triangles[i]->points[j+1]->y - yMin) * normalizeRatio * windowSize, 

                  tempBuffer, Triangles[i]->points[j]->r, Triangles[i]->points[j]->g, Triangles[i]->points[j]->b,
                             Triangles[i]->points[j+1]->r, Triangles[i]->points[j+1]->g, Triangles[i]->points[j+1]->b);
          }

        lineDDA( 
                (Triangles[i]->points[2]->x - xMin) * normalizeRatio * windowSize, 
                (Triangles[i]->points[2]->y - yMin) * normalizeRatio * windowSize, 
                (Triangles[i]->points[0]->x - xMin) * normalizeRatio * windowSize, 
                (Triangles[i]->points[0]->y - yMin) * normalizeRatio * windowSize, 

                tempBuffer, Triangles[i]->points[2]->r, Triangles[i]->points[2]->g, Triangles[i]->points[2]->b,
                           Triangles[i]->points[0]->r, Triangles[i]->points[0]->g, Triangles[i]->points[0]->b );

        for(int i = 0; i < windowSize; i++) {
          rasterize(i,tempBuffer);
        }

        for(int i = 0; i < windowSize * windowSize * 4 * 3 - 2; ++i) {
          if(tempBuffer[i] != 0) {
            mainBuffer[i] = tempBuffer[i];
          }
        }
      }
    }
  }

    void normalize_XZ() {
    float xMin = 9999; float xMax = -9999; float yMin = 9999; float yMax = -9999;
    for(int i = 0; i < Triangles.size(); i++) {
      for(int j = 0; j < Triangles[i]->points.size(); ++j) {
        xMin = min(xMin, Triangles[i]->points[j]->x);
        xMax = max(xMax, Triangles[i]->points[j]->x);
        yMin = min(yMin, Triangles[i]->points[j]->z);
        yMax = max(yMax, Triangles[i]->points[j]->z);               
      }   
    }
    float deltaX = xMax - xMin;
    float deltaY = yMax - yMin;
    normalizeRatio = normalizeSize/max(deltaY,deltaX);

    for(int i = 0; i < Triangles.size(); ++i) {
      if(Triangles[i]->normal[1] <= 0) {
        float* tempBuffer = new float [windowSize * windowSize * 3 * 4];

        for(int j = 0; j < 2; ++j){
          lineDDA( 
                  (Triangles[i]->points[j]->x - xMin) * normalizeRatio * windowSize, 
                  (Triangles[i]->points[j]->z - yMin) * normalizeRatio * windowSize + windowSize, 
                  (Triangles[i]->points[j+1]->x - xMin) * normalizeRatio * windowSize, 
                  (Triangles[i]->points[j+1]->z - yMin) * normalizeRatio * windowSize + windowSize, 

                  tempBuffer, Triangles[i]->points[j]->r, Triangles[i]->points[j]->g, Triangles[i]->points[j]->b,
                             Triangles[i]->points[j+1]->r, Triangles[i]->points[j+1]->g, Triangles[i]->points[j+1]->b);
          }

        lineDDA( 
                (Triangles[i]->points[2]->x - xMin) * normalizeRatio * windowSize, 
                (Triangles[i]->points[2]->z - yMin) * normalizeRatio * windowSize + windowSize, 
                (Triangles[i]->points[0]->x - xMin) * normalizeRatio * windowSize, 
                (Triangles[i]->points[0]->z - yMin) * normalizeRatio * windowSize + windowSize, 

                tempBuffer, Triangles[i]->points[2]->r, Triangles[i]->points[2]->g, Triangles[i]->points[2]->b,
                           Triangles[i]->points[0]->r, Triangles[i]->points[0]->g, Triangles[i]->points[0]->b );

        for(int i = windowSize-1; i < windowSize*2; i++) {
          rasterize(i,tempBuffer);
        }

        for(int i = 0; i < windowSize * windowSize * 4 * 3 - 2; ++i) {
          if(tempBuffer[i] != 0) {
            mainBuffer[i] = tempBuffer[i];
          }
        }
      }
    }
  }

    void normalize_YZ() {
    float xMin = 9999; float xMax = -9999; float yMin = 9999; float yMax = -9999;
    for(int i = 0; i < Triangles.size(); i++) {
      for(int j = 0; j < Triangles[i]->points.size(); ++j) {
        xMin = min(xMin, Triangles[i]->points[j]->y);
        xMax = max(xMax, Triangles[i]->points[j]->y);
        yMin = min(yMin, Triangles[i]->points[j]->z);
        yMax = max(yMax, Triangles[i]->points[j]->z);               
      }   
    }
    float deltaX = xMax - xMin;
    float deltaY = yMax - yMin;
    normalizeRatio = normalizeSize/max(deltaY,deltaX);

    for(int i = 0; i < Triangles.size(); ++i) {
      if(Triangles[i]->normal[0] <= 0) {
        float* tempBuffer = new float [windowSize * windowSize * 3 * 4];

        for(int j = 0; j < 2; ++j){
          lineDDA( 
                  (Triangles[i]->points[j]->y - xMin) * normalizeRatio * windowSize + windowSize, 
                  (Triangles[i]->points[j]->z - yMin) * normalizeRatio * windowSize + windowSize, 
                  (Triangles[i]->points[j+1]->y - xMin) * normalizeRatio * windowSize + windowSize, 
                  (Triangles[i]->points[j+1]->z - yMin) * normalizeRatio * windowSize + windowSize, 

                  tempBuffer, Triangles[i]->points[j]->r, Triangles[i]->points[j]->g, Triangles[i]->points[j]->b,
                             Triangles[i]->points[j+1]->r, Triangles[i]->points[j+1]->g, Triangles[i]->points[j+1]->b);
          }

        lineDDA( 
                (Triangles[i]->points[2]->y - xMin) * normalizeRatio * windowSize + windowSize, 
                (Triangles[i]->points[2]->z - yMin) * normalizeRatio * windowSize + windowSize, 
                (Triangles[i]->points[0]->y - xMin) * normalizeRatio * windowSize + windowSize, 
                (Triangles[i]->points[0]->z - yMin) * normalizeRatio * windowSize + windowSize, 

                tempBuffer, Triangles[i]->points[2]->r, Triangles[i]->points[2]->g, Triangles[i]->points[2]->b,
                           Triangles[i]->points[0]->r, Triangles[i]->points[0]->g, Triangles[i]->points[0]->b );

        for(int i = windowSize-1; i < windowSize*2; i++) {
          rasterize(i,tempBuffer);
        }

        for(int i = 0; i < windowSize * windowSize * 4 * 3 - 2; ++i) {
          if(tempBuffer[i] != 0) {
            mainBuffer[i] = tempBuffer[i];
          }
        }
      }
    }
  }

  void normalizeAll(){
    //sortSurfaces(1);
    normalize_XY();
    //sortSurfaces(3);
    //normalize_XZ();
    //sortSurfaces(2);
    //normalize_YZ();
  }

  void clearBuffer() {
      mainBuffer = new float[windowSize * windowSize * 3 * 4];
  }

  void WomboCombo() {
    clearBuffer();
    configurationDisplayed();
    allObjectsRecompute();
    transferTriangles();
    allTrianglesRecompute();
    allVertexRecompute();
    normalizeAll();
  }
};

Scene manager(0.95);

void lineDDA (int x0, int y0, int xEnd, int yEnd, float* board, float r0, float g0, float b0, float rEnd, float gEnd, float bEnd) {

  float deltaR = rEnd - r0;
  float deltaG = gEnd - g0;
  float deltaB = bEnd - b0;

  int dx = xEnd - x0, dy = yEnd - y0, steps, k;
  float xIncrement, yIncrement, x = x0, y = y0;

  if(fabs(dx) > fabs(dy)){
    steps = fabs(dx);
  } else {
    steps = fabs(dy);
  }
  xIncrement = float(dx) / float(steps);
  yIncrement = float(dy) / float(steps);
  makePix(round(x), round(y), board, r0 , g0, b0);

  for(int k = 0; k < steps; ++k){
    x += xIncrement;
    y += yIncrement;
    makePix(round(x),round(y), board, r0 + (k+1)*(deltaR/steps), g0 + (k+1)*(deltaG/steps) , b0 + (k+1)*(deltaB/steps));
  }
}

void rasterize(int y, float* board) {  
  int x0,xEnd;
  x0 = 0;
  xEnd = windowSize * 2 - 1;
  float r0,g0,b0,rEnd,gEnd,bEnd;
  int state = 0;
  int temp;
  for(int i = 0; i < windowSize * 2 ; i++) {
    temp = 3 * windowSize * 2 * y + i * 3;
    switch(state) {
      case 0 : 
        if(board[temp] != 0) {
          x0 = i;
          r0 = board[temp];
          g0 = board[temp+1];
          b0 = board[temp+2];
          state = 1;
        }
      break;
      case 1 :
        if(board[temp] == 0) {
          state = 2;
        }
      break;
      case 2:
        if(board[temp] != 0) {
          xEnd = i;
          rEnd = board[temp];
          gEnd = board[temp+1];
          bEnd = board[temp+2];
          state = 3;
        }
      break;
      case 3: // do nothing
      break;
    } 
  }

  float deltaR = rEnd-r0;
  float deltaG = gEnd-g0;
  float deltaB = bEnd-b0; 

  if(state == 3){
    for(int i = x0 + 1; i < xEnd; ++i) {
      makePix(i, y, board, r0 + (i-x0) * deltaR/float(xEnd - x0), g0 + (i-x0) * deltaG/float(xEnd - x0), b0 + (i-x0) * deltaB/float(xEnd - x0));
    }
  }
}

void makePix(int x, int y, float* board, float r, float g, float b) {
  board[3 * windowSize * y * 2 + x * 3    ] = r;
  board[3 * windowSize * y * 2 + x * 3 + 1] = g;
  board[3 * windowSize * y * 2 + x * 3 + 2] = b;
}

void init() {
  glClearColor(0.0,0.0,0.0,0.0);
}

void display(void) 
{
  glClear(GL_COLOR_BUFFER_BIT); 
  glLoadIdentity(); 
  glColor3f(1.0,0.0,1.0);
  glDrawPixels(windowSize * 2, windowSize * 2, GL_RGB, GL_FLOAT, mainBuffer);
  glFlush(); 
}

void keyboard(unsigned char c, int x, int y) {

  float speedfactor = 0.17;
  switch(c) {
    case '1' : objectNum = 0;  break;
    case '2' : objectNum = 1;  break;
    case '3' : objectNum = 2;  break;

    case 'a' :
      Objects[objectNum].translation(-0.3*speedfactor,0,0);
      manager.WomboCombo();
      glutPostRedisplay();
      break;

    case 'd' :
      Objects[objectNum].translation(0.3*speedfactor,0,0);
      manager.WomboCombo();
      glutPostRedisplay();
      break;

    case 'w' :
      Objects[objectNum].translation(0,0.3*speedfactor,0);
      manager.WomboCombo();
      glutPostRedisplay();
      break;

    case 's' :
      Objects[objectNum].translation(0,-0.3*speedfactor,0);
      manager.WomboCombo();
      glutPostRedisplay();
      break;

    case 'z' :
      Objects[objectNum].scale(1 - 0.5 * speedfactor);
      manager.WomboCombo();
      glutPostRedisplay();
      break;

    case 'c' :
      Objects[objectNum].scale(1 + 0.5 * speedfactor);
      manager.WomboCombo();
      glutPostRedisplay();
      break;

    case 'i' :
      Objects[objectNum].rotationUpInto(-10*speedfactor,1);
      manager.WomboCombo();
      glutPostRedisplay();      
      break;  

    case 'k' :
      Objects[objectNum].rotationUpInto(10*speedfactor,1);
      manager.WomboCombo();
      glutPostRedisplay();      
      break; 

    case 'j' :
      Objects[objectNum].rotationUpInto(-10*speedfactor,3);
      manager.WomboCombo();
      glutPostRedisplay();      
      break;  

    case 'l' :
      Objects[objectNum].rotationUpInto(10*speedfactor,3);
      manager.WomboCombo();
      glutPostRedisplay();      
      break;

    case 'p': 
      ofstream output;
      output.open("vertexNormals.txt");
      output << Objects.size() << endl;
      for(int i = 0; i < Objects.size(); ++i) {
        Objects[i].computeVertexNormal();
        output << "Object #" << i << endl;
        for(int j = 0; j < Objects[i].points.size(); ++j) {
          output << "Vertex #" << j << endl;
          output << Objects[i].points[j].normal[0] << " " << Objects[i].points[j].normal[1] << " " << Objects[i].points[j].normal[2] << endl;
        }
      }
      output.close();
      break;
  }
}


int main(int argc, char* argv[]) {
  windowSize = 350;

  cout << "Welcome to Proj3" << endl;
  int option;

  cout << "Please select: \n1) Load default perimeter values \n2) Define perimeter values. \n";
  cin >> option;
  if(option == 2) {
    cout << "enter Ka: ";
    cin >> Ka;
    cout << "\nenter Kd: ";
    cin >> Kd;
    cout << "\nenter Ks: ";
    cin >> Ks;
    cout << "\nenter n: ";
    cin >> n;
    cout << "\nenter IA: ";
    cin >> IA;
    cout << "\nenter IL: ";
    cin >> IL;
  } else {
    IL = 0.6;
    Kd = 0.9;
    Ks = 0.4;
    n = 3;
    IA = 0.3;
    IL = 0.4;
  }

  lightSource[0] = 0;
  lightSource[1] = 0.1;
  lightSource[2] = 0.05;

  mainBuffer = new float[windowSize * windowSize * 4 * 3];

  string line;
  ifstream myfile;
  myfile.open("Shapes.txt");

  int numObj, numSurf, numPts;
  float x,y,z;
  int p1,p2,p3;
  getline(myfile,line);

  numObj = atoi(line.c_str());
  while(numObj > 0){
    Object temp;
    getline(myfile,line);
    numPts = atoi(line.c_str());
    while(numPts > 0){
      getline(myfile,line);
      stringstream(line) >> x >> y >> z;
      temp.addVertex(x,y,z);
      numPts--;
    }

    getline(myfile,line);
    numSurf = atoi(line.c_str());

    while(numSurf > 0){
      getline(myfile,line);
      stringstream haha(line);
      haha >> p1 >> p2 >> p3;
      temp.addSurface(p1,p2,p3);
      numSurf--;
    }

    temp.computeValues();
    Objects.push_back(temp);
    numObj--;
  }


  for(int i = 0; i < Objects.size(); ++i) {
    Objects[i].generateSurfaces();
    Objects[i].computeVertexNormal();
    Objects[i].computeVertexIntensity();

  }

  myfile.close();
  manager.allObjectsRecompute();
  manager.transferTriangles();
  manager.allTrianglesRecompute();
  manager.normalizeAll();

  glutInit(&argc, argv);           
  glutInitDisplayMode(GLUT_SINGLE);
  glutInitWindowSize (windowSize * 2, windowSize * 2);  
  glutInitWindowPosition (100,100);
  glutCreateWindow ("Project 3 ");
  init();
  glutKeyboardFunc(keyboard);
  glutDisplayFunc(display);
  glutMainLoop();
}
