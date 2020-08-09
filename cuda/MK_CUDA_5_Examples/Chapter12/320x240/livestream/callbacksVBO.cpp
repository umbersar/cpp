#include <GL/glew.h>
#include <GL/freeglut.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include <helper_timer.h>
#include <helper_cuda.h>
#include <helper_cuda_gl.h>

// The user must create the following routines:
void initCuda(int argc, char** argv);
void runCuda();
void renderCuda(int);

// Callback variables
extern float animTime;
extern int sleepTime, sleepInc;
int drawMode=GL_TRIANGLE_FAN; // the default draw mode
int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;
float rotate_x = 0.0, rotate_y = 0.0;
float translate_z = -3.0;

// some initial values for face segmentation
int PureG[2]={62,89}, PureR[2]={112,145};
int doWave=1, doSkin=0, doSobel=0;

// GLUT callbacks display, keyboard, mouse
void display()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // set view matrix
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(0.0, 0.0, translate_z);
  glRotatef(rotate_x, 1.0, 0.0, 0.0);
  glRotatef(rotate_y, 0.0, 1.0, 0.0);

  runCuda(); // run CUDA kernel to generate vertex positions

  renderCuda(drawMode); // render the data
  
  glutSwapBuffers();
  glutPostRedisplay();
  
  // slow the rendering when the GPU is too fast.
  if(sleepTime) usleep(sleepTime);
  animTime += 0.01;
}

void keyboard(unsigned char key, int x, int y)
{
  switch(key) {
  case('q') : case(27) : // exit
    exit(0);
    break;
  case 'd': case 'D': // Drawmode
    switch(drawMode) {
    case GL_POINTS: drawMode = GL_LINE_STRIP; break;
    case GL_LINE_STRIP: drawMode = GL_TRIANGLE_FAN; break;
    default: drawMode=GL_POINTS;
    } break;
  case 'S': // Slow the simulation down
    sleepTime += sleepInc;
    break;
  case 's': // Speed the simulation up
    sleepTime = (sleepTime > 0)?sleepTime -= sleepInc:0;
    break;
  case 'z': doWave = (doWave > 0)?0:1; break;
  case 'x': doSkin = (doSkin > 0)?0:1; break;
  case 'c': doSobel = (doSobel > 0)?0:1; break;
  case 'R': PureR[1]++; if(PureR[1] > 255) PureR[1]=255; break;
  case 'r': PureR[1]--; if(PureR[1] <= PureR[0]) PureR[1]++; break;
  case 'E': PureR[0]++; if(PureR[0] >= PureR[1]) PureR[0]--; break;
  case 'e': PureR[0]--; if(PureR[0] <= 0 ) PureR[0]=0; break;
  case 'G': PureG[1]++; if(PureG[1] > 255) PureG[1]=255; break;
  case 'g': PureG[1]--; if(PureG[1] <= PureG[0]) PureG[1]++; break;
  case 'F': PureG[0]++; if(PureG[0] >= PureG[1]) PureG[0]--; break;
  case 'f': PureG[0]--; if(PureG[0] <= 0 ) PureG[0]=0; break;
  }
  fprintf(stderr,"PureG[0] %d PureG[1] %d PureR[0] %d PureR[1] %d\n",
	  PureG[0],PureG[1],PureR[0],PureR[1]);
  glutPostRedisplay();
}

void mouse(int button, int state, int x, int y)
{
  if (state == GLUT_DOWN) {
    mouse_buttons |= 1<<button;
  } else if (state == GLUT_UP) {
    mouse_buttons = 0;
  }
  
  mouse_old_x = x;
  mouse_old_y = y;
  glutPostRedisplay();
}

void motion(int x, int y)
{
  float dx, dy;
  dx = x - mouse_old_x;
  dy = y - mouse_old_y;
  
  if (mouse_buttons & 1) {
    rotate_x += dy * 0.2;
    rotate_y += dx * 0.2;
  } else if (mouse_buttons & 4) {
    translate_z += dy * 0.01;
  }
 rotate_x = (rotate_x < -60.)?-60.:(rotate_x > 60.)?60:rotate_x;
 rotate_y = (rotate_y < -60.)?-60.:(rotate_y > 60.)?60:rotate_y;
  
  mouse_old_x = x;
  mouse_old_y = y;
}


