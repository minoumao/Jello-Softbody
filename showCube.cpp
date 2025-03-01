/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "showCube.h"
#include <vector>

int pointMap(int side, int i, int j)
{
  int r;

  switch (side)
  {
  case 1: //[i][j][0] bottom face
    r = 64 * i + 8 * j;
    break;
  case 6: //[i][j][7] top face
    r = 64 * i + 8 * j + 7;
    break;
  case 2: //[i][0][j] front face
    r = 64 * i + j;
    break;
  case 5: //[i][7][j] back face
    r = 64 * i + 56 + j;
    break;
  case 3: //[0][i][j] left face
    r = 8 * i + j;
    break;
  case 4: //[7][i][j] right face
    r = 448 + 8 * i + j;
    break;
  }

  return r;
}

void showCube(struct world * jello)
{
  int i,j,k,ip,jp,kp;
  point r1,r2,r3; // aux variables
  
  /* normals buffer and counter for Gourad shading*/
  struct point normal[8][8];
  int counter[8][8];

  int face;
  double faceFactor, length;

  if (fabs(jello->p[0][0][0].x) > 10)
  {
    printf ("Your cube somehow escaped way out of the box.\n");
    exit(0);
  }

  
  #define NODE(face,i,j) (*((struct point * )(jello->p) + pointMap((face),(i),(j))))

  
  #define PROCESS_NEIGHBOUR(di,dj,dk) \
    ip=i+(di);\
    jp=j+(dj);\
    kp=k+(dk);\
    if\
    (!( (ip>7) || (ip<0) ||\
      (jp>7) || (jp<0) ||\
    (kp>7) || (kp<0) ) && ((i==0) || (i==7) || (j==0) || (j==7) || (k==0) || (k==7))\
       && ((ip==0) || (ip==7) || (jp==0) || (jp==7) || (kp==0) || (kp==7))) \
    {\
      glVertex3f(jello->p[i][j][k].x,jello->p[i][j][k].y,jello->p[i][j][k].z);\
      glVertex3f(jello->p[ip][jp][kp].x,jello->p[ip][jp][kp].y,jello->p[ip][jp][kp].z);\
    }\

 
  if (viewingMode==0) // render wireframe
  {
    glLineWidth(1);
    glPointSize(5);
    glDisable(GL_LIGHTING);
    for (i=0; i<=7; i++)
      for (j=0; j<=7; j++)
        for (k=0; k<=7; k++)
        {
          if (i*j*k*(7-i)*(7-j)*(7-k) != 0) // not surface point
            continue;

          glBegin(GL_POINTS); // draw point
            glColor4f(0,0,0,0);  
            glVertex3f(jello->p[i][j][k].x,jello->p[i][j][k].y,jello->p[i][j][k].z);        
          glEnd();

          //
          //if ((i!=7) || (j!=7) || (k!=7))
          //  continue;

          glBegin(GL_LINES);      
          // structural
          if (structural == 1)
          {
            glColor4f(0,0,1,1);
            PROCESS_NEIGHBOUR(1,0,0);
            PROCESS_NEIGHBOUR(0,1,0);
            PROCESS_NEIGHBOUR(0,0,1);
            PROCESS_NEIGHBOUR(-1,0,0);
            PROCESS_NEIGHBOUR(0,-1,0);
            PROCESS_NEIGHBOUR(0,0,-1);
          }
          
          // shear
          if (shear == 1)
          {
            glColor4f(0,1,0,1);
            PROCESS_NEIGHBOUR(1,1,0);
            PROCESS_NEIGHBOUR(-1,1,0);
            PROCESS_NEIGHBOUR(-1,-1,0);
            PROCESS_NEIGHBOUR(1,-1,0);
            PROCESS_NEIGHBOUR(0,1,1);
            PROCESS_NEIGHBOUR(0,-1,1);
            PROCESS_NEIGHBOUR(0,-1,-1);
            PROCESS_NEIGHBOUR(0,1,-1);
            PROCESS_NEIGHBOUR(1,0,1);
            PROCESS_NEIGHBOUR(-1,0,1);
            PROCESS_NEIGHBOUR(-1,0,-1);
            PROCESS_NEIGHBOUR(1,0,-1);

            PROCESS_NEIGHBOUR(1,1,1)
            PROCESS_NEIGHBOUR(-1,1,1)
            PROCESS_NEIGHBOUR(-1,-1,1)
            PROCESS_NEIGHBOUR(1,-1,1)
            PROCESS_NEIGHBOUR(1,1,-1)
            PROCESS_NEIGHBOUR(-1,1,-1)
            PROCESS_NEIGHBOUR(-1,-1,-1)
            PROCESS_NEIGHBOUR(1,-1,-1)
          }
          
          // bend
          if (bend == 1)
          {
            glColor4f(1,0,0,1);
            PROCESS_NEIGHBOUR(2,0,0);
            PROCESS_NEIGHBOUR(0,2,0);
            PROCESS_NEIGHBOUR(0,0,2);
            PROCESS_NEIGHBOUR(-2,0,0);
            PROCESS_NEIGHBOUR(0,-2,0);
            PROCESS_NEIGHBOUR(0,0,-2);
          }           
          glEnd();
        }
    glEnable(GL_LIGHTING);
  }
  
  else
  {
    glPolygonMode(GL_FRONT, GL_FILL); 
    
    for (face=1; face <= 6; face++) 
      // face == face of a cube
      // 1 = bottom, 2 = front, 3 = left, 4 = right, 5 = far, 6 = top
    {
      
      if ((face==1) || (face==3) || (face==5))
        faceFactor=-1; // flip orientation
      else
        faceFactor=1;
      

      for (i=0; i <= 7; i++) // reset buffers
        for (j=0; j <= 7; j++)
        {
          normal[i][j].x=0;normal[i][j].y=0;normal[i][j].z=0;
          counter[i][j]=0;
        }

      /* process triangles, accumulate normals for Gourad shading */
  
      for (i=0; i <= 6; i++)
        for (j=0; j <= 6; j++) // process block (i,j)
        {
          pDIFFERENCE(NODE(face,i+1,j),NODE(face,i,j),r1); // first triangle
          pDIFFERENCE(NODE(face,i,j+1),NODE(face,i,j),r2);
          CROSSPRODUCTp(r1,r2,r3); pMULTIPLY(r3,faceFactor,r3);
          pNORMALIZE(r3);
          pSUM(normal[i+1][j],r3,normal[i+1][j]);
          counter[i+1][j]++;
          pSUM(normal[i][j+1],r3,normal[i][j+1]);
          counter[i][j+1]++;
          pSUM(normal[i][j],r3,normal[i][j]);
          counter[i][j]++;

          pDIFFERENCE(NODE(face,i,j+1),NODE(face,i+1,j+1),r1); // second triangle
          pDIFFERENCE(NODE(face,i+1,j),NODE(face,i+1,j+1),r2);
          CROSSPRODUCTp(r1,r2,r3); pMULTIPLY(r3,faceFactor,r3);
          pNORMALIZE(r3);
          pSUM(normal[i+1][j],r3,normal[i+1][j]);
          counter[i+1][j]++;
          pSUM(normal[i][j+1],r3,normal[i][j+1]);
          counter[i][j+1]++;
          pSUM(normal[i+1][j+1],r3,normal[i+1][j+1]);
          counter[i+1][j+1]++;
        }

      
        /* the actual rendering */
        for (j=1; j<=7; j++) 
        {

          if (faceFactor  > 0)
            glFrontFace(GL_CCW); // the usual definition of front face
          else
            glFrontFace(GL_CW); // flip definition of orientation
         
          glBegin(GL_TRIANGLE_STRIP);
          for (i=0; i<=7; i++)
          {
            glNormal3f(normal[i][j].x / counter[i][j],normal[i][j].y / counter[i][j],
              normal[i][j].z / counter[i][j]);
            glVertex3f(NODE(face,i,j).x, NODE(face,i,j).y, NODE(face,i,j).z);
            glNormal3f(normal[i][j-1].x / counter[i][j-1],normal[i][j-1].y/ counter[i][j-1],
              normal[i][j-1].z / counter[i][j-1]);
            glVertex3f(NODE(face,i,j-1).x, NODE(face,i,j-1).y, NODE(face,i,j-1).z);
          }
          glEnd();
        }
        
        
    }  
  } // end for loop over faces
  glFrontFace(GL_CCW);
}

void showBoundingBox()
{
  int i,j;

  glColor4f(0.6,0.6,0.6,0);

  glBegin(GL_LINES);

  // front face
  for(i=-2; i<=2; i++)
  {
    glVertex3f(i,-2,-2);
    glVertex3f(i,-2,2);
  }
  for(j=-2; j<=2; j++)
  {
    glVertex3f(-2,-2,j);
    glVertex3f(2,-2,j);
  }

  // back face
  for(i=-2; i<=2; i++)
  {
    glVertex3f(i,2,-2);
    glVertex3f(i,2,2);
  }
  for(j=-2; j<=2; j++)
  {
    glVertex3f(-2,2,j);
    glVertex3f(2,2,j);
  }

  // left face
  for(i=-2; i<=2; i++)
  {
    glVertex3f(-2,i,-2);
    glVertex3f(-2,i,2);
  }
  for(j=-2; j<=2; j++)
  {
    glVertex3f(-2,-2,j);
    glVertex3f(-2,2,j);
  }

  // right face
  for(i=-2; i<=2; i++)
  {
    glVertex3f(2,i,-2);
    glVertex3f(2,i,2);
  }
  for(j=-2; j<=2; j++)
  {
    glVertex3f(2,-2,j);
    glVertex3f(2,2,j);
  }
  
  glEnd();

  return;
}

void showGroundPlane()
{
    int i = -1;
    int j = -2;

    // first stagger
    while (i < 2)
    {
        while (j < 2)
        {
            // set tile color based on modulus
            if (j % 2 == 0)
                glColor4f(0.1, 0.1, 0.1, 1.0);
            else
                glColor4f(0.0, 0.0, 0.0, 0.0);

            // draw two squares at the same time. Triangles CCW
            glBegin(GL_TRIANGLES);
            // triangle 1
            glVertex3f(i, j, -2);
            glVertex3f(i - 1, j + 1, -2);
            glVertex3f(i - 1, j, -2);

            // triangle 2
            glVertex3f(i, j, -2);
            glVertex3f(i, j + 1, -2);
            glVertex3f(i - 1, j + 1, -2);
            glEnd();

            glBegin(GL_TRIANGLES);
            // triangle 1
            glVertex3f(i + 2, j, -2);
            glVertex3f(i + 1, j + 1, -2);
            glVertex3f(i + 1, j, -2);

            // triangle 2
            glVertex3f(i + 2, j, -2);
            glVertex3f(i + 2, j + 1, -2);
            glVertex3f(i + 1, j + 1, -2);
            glEnd();
            j++;
        }
        i++;
    }

    // second stagger
    i = -1;
    j = -2;
    while (i < 2)
    {
        while (j < 2)
        {
            // set tile color based on modulus
            if (j % 2 == 0)
                glColor4f(0.0, 0.0, 0.0, 0.0);
            else
                glColor4f(0.1, 0.1, 0.1, 1.0);

            // draw two squares at the same time. Triangles CCW
            glBegin(GL_TRIANGLES);
            // triangle 1
            glVertex3f(i + 1, j, -2);
            glVertex3f(i, j + 1, -2);
            glVertex3f(i, j, -2);

            // triangle 2
            glVertex3f(i + 1, j, -2);
            glVertex3f(i + 1, j + 1, -2);
            glVertex3f(i, j + 1, -2);
            glEnd();

            glBegin(GL_TRIANGLES);
            // triangle 1
            glVertex3f(i + 3, j, -2);
            glVertex3f(i + 2, j + 1, -2);
            glVertex3f(i + 2, j, -2);

            // triangle 2
            glVertex3f(i + 3, j, -2);
            glVertex3f(i + 3, j + 1, -2);
            glVertex3f(i + 2, j + 1, -2);
            glEnd();
            j++;
        }
        i++;
    }

}


void showInclinePlane(struct world* jello)
{
    // define 8 points
    // front surface
    double _frontLeft;
    double _frontRight;

    // back surface
    double _backLeft;
    double _backRight;

    // up surface
    double _upFront;
    double _upBack;

    // bottom surface
    double _bottomFront;
    double _bottomBack;

    // point storage
    point _points[8];
    bool _intersectionPoints[8];

    // Initalize
    for (int i = 0; i < 8; i++)
    {
        _intersectionPoints[i] = false;
    }

    // front left
    _frontLeft = -((jello->a * (-2.0) + jello->c * (-2.0) + jello->d) / jello->b);
    _points[0].x = -2;
    _points[0].z = -2;
    if ((_frontLeft >= -2 && _frontLeft <= 2))
    {
        _intersectionPoints[0] = true;
        _points[0].y = _frontLeft;
    }

    // front right
    _frontRight = -((jello->a * (2.0) + jello->c * (-2.0) + jello->d) / jello->b);
    _points[1].x = 2;
    _points[1].z = -2;
    if ((_frontRight >= -2 && _frontRight <= 2))
    {
        _intersectionPoints[1] = true;
        _points[1].y = _frontRight;
    }

    // back left
    _backLeft = -((jello->a * (-2.0) + jello->c * (2.0) + jello->d) / jello->b);
    _points[2].x = -2;
    _points[2].z = 2;
    if ((_backLeft >= -2 && _backLeft <= 2))
    {
        _intersectionPoints[2] = true;
        _points[2].y = _backLeft;
    }

    // back right
    _backRight = -((jello->a * (2.0) + jello->c * (2.0) + jello->d) / jello->b);
    _points[3].x = 2;
    _points[3].z = 2;
    if ((_backRight >= -2 && _backRight <= 2))
    {
        _intersectionPoints[3] = true;
        _points[3].y = _backRight;
    }

    // up front
    _upFront = -((jello->b * (2.0) + jello->c * (-2.0) + jello->d) / jello->a);
    _points[4].y = 2;
    _points[4].z = -2;
    if ((_upFront > -2 && _upFront < 2))
    {
        _intersectionPoints[4] = true;
        _points[4].x = _upFront;
    }

    // up back
    _upBack = -((jello->b * (2.0) + jello->c * (2.0) + jello->d) / jello->a);
    _points[5].y = 2;
    _points[5].z = 2;
    if ((_upBack > -2 && _upBack < 2))
    {
        _intersectionPoints[5] = true;
        _points[5].x = _upBack;
    }

    // bottom front
    _bottomFront = -((jello->b * (-2.0) + jello->c * (-2.0) + jello->d) / jello->a);
    _points[6].y = -2;
    _points[6].z = -2;
    if ((_bottomFront > -2 && _bottomFront < 2))
    {
        _intersectionPoints[6] = true;
        _points[6].x = _bottomFront;
    }

    // bottom back
    _bottomBack = -((jello->b * (-2.0) + jello->c * (2.0) + jello->d) / jello->a);
    _points[7].y = -2;
    _points[7].z = 2;
    if ((_bottomBack > -2 && _bottomBack < 2))
    {
        _intersectionPoints[7] = true;
        _points[7].x = _bottomBack;
    }

    // gather points
    int pointCount = 0;
    std::vector<point> drawPoint;
    for (int i = 0; i < 8; i++)
    {
        if (_intersectionPoints[i])
        {
            pointCount++;
            drawPoint.push_back(_points[i]);
        }
    }

    // draw triangle
    if (pointCount == 3)
    {
        glBegin(GL_POLYGON);

        glColor4f(0.2, 0.2, 0.2, 1);
        glVertex3f(drawPoint[2].x, drawPoint[2].y, drawPoint[2].z);

        glColor4f(0.2, 0.2, 0.2, 1);
        glVertex3f(drawPoint[1].x, drawPoint[1].y, drawPoint[1].z);

        glColor4f(0.2, 0.2, 0.2, 1);
        glVertex3f(drawPoint[0].x, drawPoint[0].y, drawPoint[0].z);

        glEnd();
    }

    // draw quad
    if(pointCount == 4)
    {
        glBegin(GL_QUADS);

        glColor4f(0.2, 0.2, 0.2, 1);
        glVertex3f(_points[3].x, _backRight, _points[3].z);
        glVertex3f(_points[0].x, _frontLeft, _points[0].z);

        glColor4f(0.2, 0.2, 0.2, 1);
        glVertex3f(_points[1].x, _frontRight, _points[1].z);
        glVertex3f(_points[2].x, _backLeft, _points[2].z);

        glEnd();
    }
}