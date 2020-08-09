// simpleGLmain (Rob Farber)
// modified to build under CUDA-5
#include <stdio.h>

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include <helper_timer.h>

// GLUT specific contants
const unsigned int window_width = 512;
const unsigned int window_height = 512;

// The user must create the following routines:
void initCuda(int argc, char** argv);
bool initGL(int argc, char** argv);
void fpsDisplay(), display();
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);

StopWatchInterface *timer = NULL; // a timer for FPS calculations
int sleepTime=0, sleepInc=100;

// Main program
int main(int argc, char** argv)
{
  // Create the timer
  sdkCreateTimer(&timer);
  
  if (false == initGL(argc, argv)) { return false; }


  initCuda(argc, argv);

  // register callbacks
  glutDisplayFunc(fpsDisplay);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  
  // start rendering mainloop
  glutMainLoop();
  
  // clean up
  cudaThreadExit();
  cudaDeviceReset();
}

// Simple method to display the Frames Per Second in the window title
void computeFPS()
{
  static int fpsCount=0;
  static int fpsLimit=100;

  fpsCount++;
  
  if (fpsCount == fpsLimit) {
    char fps[256];
    float ifps = 1.f / (sdkGetAverageTimerValue(&timer) / 1000.f);
    if(sleepTime)
      sprintf(fps, "CUDA Interop (Rob Farber): %3.1f fps sleepTime %3.1f ms ",
	      ifps, sleepTime/1000.);  
    else
      sprintf(fps, "CUDA Interop (Rob Farber): %3.1f fps ", ifps);  
    
    glutSetWindowTitle(fps);
    fpsCount = 0; 
    
    sdkResetTimer(&timer);  
  }
}

void fpsDisplay()
{
  sdkStartTimer(&timer);  
  display();
  sdkStopTimer(&timer);
  computeFPS();
}

float animTime = 0.0;    // time the animation has been running

// Initialize OpenGL window
bool initGL(int argc, char **argv)
{
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
  glutInitWindowSize(window_width, window_height);
  glutCreateWindow("Cuda GL Interop Demo (adapted from NVIDIA's simpleGL");
  glutDisplayFunc(fpsDisplay);
  glutKeyboardFunc(keyboard);
  glutMotionFunc(motion);
  
  // initialize necessary OpenGL extensions
  glewInit();
  if (! glewIsSupported("GL_VERSION_2_0 ")) {
    fprintf(stderr, "ERROR: Support for necessary OpenGL extensions missing.");
    fflush(stderr);
    return false;
  }
  
  // default initialization
  glClearColor(0.0, 0.0, 0.0, 1.0);
  glDisable(GL_DEPTH_TEST);
  
  // viewport
  glViewport(0, 0, window_width, window_height);
  
  // set view matrix
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // projection
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0, (GLfloat)window_width/(GLfloat) window_height,0.10, 10.0);
  
  return true;
}

