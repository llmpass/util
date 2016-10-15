#ifndef _SCRSHOT_H
#define _SCRSHOT_H

#include <iostream>

void printToPPM(const int width, const int height, const char* sn) {
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glFlush();
  glutSwapBuffers();
  char name[1024];
  strcat(strcpy(name,sn),".ppm");
  ofstream ofs(name);
  unsigned char *image = new unsigned char[width*height*3];
  glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,image);
  // vertical axes is upside-down
  unsigned char *imageU = new unsigned char[width*height*3];
  for (int i=0; i<height; ++i)
    for (int j=0; j<width; ++j) 
      for (int k=0; k<3; ++k)  
        imageU[i*width*3+j*3+k] = 
          image[(height-1-i)*width*3+j*3+k];
  ofs << "P6" << endl << width << " " << height << endl << 255 << endl;
  ofs.write((char*)imageU,width*height*3);
  ofs.close();
}

#endif
