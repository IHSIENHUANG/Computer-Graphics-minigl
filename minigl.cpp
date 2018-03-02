/*gl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>
#include <limits.h>
using namespace std;
#define DEBUG 0
#define PI 3.141592
/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h
vec3 RGB;
struct Vertex
{
    vec4 position;
    vec3 color;
    Vertex(vec4 pos,vec3 co)
    {
        position=pos;
        color=co;    
    };
    Vertex(){};
};
struct Triangle
{
  Vertex A,B,C;
  Triangle(Vertex a, Vertex b, Vertex c)
  {
    A=a;
    B=b;
    C=c;          
  }; 
  Triangle()
  {
  };
};
vector<Vertex> list_vertex;
vector<Triangle> list_triangle;  
MGLpoly_mode type;
mat4 matrix_4 = {{1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1}};
mat4 identity_matrix =  {{1.0f,0,0,0,0,1.0f,0,0,0,0,1.0f,0,0,0,0,1.0f}};
MGLmatrix_mode pro_mode;
vector<mat4> Projection_stack(1,identity_matrix);
vector<mat4> Model_stack(1,identity_matrix);
vector<int> z_buffer(500*500,INT_MAX);
/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}


/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void Rasterize_Triangle(const Triangle& tri,int width,int height,MGLpixel* data)
{
  // cout<<"width "<<width<<" height "<<height<<endl; 
   for(int i = 0 ;i<width;i++)
    {
        for(int j = 0 ; j < height ;j++)
        {
                     
                double fi,fj;
                //fj = (j+1)*height/2.0-0.5;
                fi=i;
                fj=j;
                
                Vertex T_A =  tri.A;
                Vertex T_B =  tri.B;
                Vertex T_C =  tri.C;
                //we need to count w y z by devide by w 
                double a_x = T_A.position[0]/T_A.position[3],a_y = T_A.position[1]/T_A.position[3];
                double b_x = T_B.position[0]/T_B.position[3],b_y = T_B.position[1]/T_B.position[3];
                double c_x = T_C.position[0]/T_C.position[3],c_y = T_C.position[1]/T_C.position[3];
                    
                a_x=(a_x+1)*width/2.0-0.5;
                b_x=(b_x+1)*width/2.0-0.5;
                c_x=(c_x+1)*width/2.0-0.5;
                a_y=(a_y+1)*height/2.0-0.5;
                b_y=(b_y+1)*height/2.0-0.5;
                c_y=(c_y+1)*height/2.0-0.5;
               // cout<<a_x <<" "<<a_y<<" "<<b_x<<" "<<b_y<<" "<<c_x<<" "<<c_y<<endl;
                double area_ABC =a_x*(b_y-c_y)+a_y*(c_x-b_x)+(b_x*c_y-b_y*c_x);

                double area_PBC =fi*(b_y-c_y)+fj*(c_x-b_x)+(b_x*c_y-b_y*c_x);
                double area_APC =a_x*(fj-c_y)+a_y*(c_x-fi)+(fi*c_y-fj*c_x);
                double area_ABP =a_x*(b_y-fj)+a_y*(fi-b_x)+(b_x*fj-b_y*fi);
                double arfa,beta,gama;
                arfa = area_PBC/area_ABC;
                beta = area_APC/area_ABC;
                gama = area_ABP/area_ABC;
                double z = arfa*T_A.position[2]+beta*T_B.position[2]+gama*T_C.position[2];//stroe z position
                int flag=0;
                if(arfa>=0 && arfa <=1) flag++;
                if(beta>=0 && beta <=1) flag++;
                if(gama>=0 && gama <=1) flag++;
                if(flag==3)
                {
                    //data[i+width*j] = Make_Pixel(255,255,255);
                   // cout<<"arfe beta gama "<<arfa+beta+gama<<endl;
                    if(data[i+width*j]==255 || z<z_buffer[i+width*j] )//maintain a global z_buffer
                    {
                        vec3 bay_color;
                        for(int k = 0 ; k <3;k++)
                        {
                            bay_color[k]=arfa*T_A.color[k]+beta*T_B.color[k]+gama*T_C.color[k];    
                        }
                        
                        
                        data[i+width*j] = Make_Pixel(bay_color[0]*255.0f,bay_color[1]*255.0f,bay_color[2]*255.0f);
                        z_buffer[i+width*j] = z;
                        
                    }
                }
        }
    }
}


void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{

   for(unsigned int k = 0 ; k < list_triangle.size();k++)//calculate each element in dif tri
   {
                Rasterize_Triangle(list_triangle[k],width,height,data);                   
   }     
    list_triangle.clear();


}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
    type =mode;
}


/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
    if(type == MGL_TRIANGLES)
    {
        for(unsigned int i = 0 ; i< list_vertex.size()/3*3;i+=3)
        {
           Triangle tri;
           tri.A = list_vertex[i];
           tri.B = list_vertex[i+1];
           tri.C = list_vertex[i+2];
           list_triangle.push_back(tri); 
        }
    }
    else if(type == MGL_QUADS)//each QUAD contains two triangle;
    {
        for(unsigned int i =0;i<list_vertex.size()/4*4;i+=4)
        {
            Triangle tri1;
            Triangle tri2;
                
            tri1.A=list_vertex[i];
            tri1.B=list_vertex[i+1];
            tri1.C=list_vertex[i+2];
            //-----
            tri2.A=list_vertex[i];
            tri2.B=list_vertex[i+2];
            tri2.C=list_vertex[i+3];
                
            list_triangle.push_back(tri1);
            list_triangle.push_back(tri2);
        }
    }
    list_vertex.clear();//clean the list of vertices 
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
    mglVertex3(x,y,0);
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
    vec4 pos;
    pos[0]=x;
    pos[1]=y;
    pos[2]=z;
    pos[3]=1.0;
    if(pro_mode==MGL_PROJECTION  )
    {
        pos = Projection_stack.back() * pos;
#if DEBUG        
        cout<<"this is projection mode"<<endl;
#endif     
    }
    else if(pro_mode==MGL_MODELVIEW)
    {
        pos = Model_stack.back()*pos;    
    }
    Vertex vertex(pos,RGB);
    list_vertex.push_back(vertex);
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{

    // Projection_stack.clear();         
    // if(!Projection_stack.empty())
//	printf("this is not clena\n");  
  //   Projection_stack.push_back(identity_matrix);
    // Model_stack.clear();
    // Model_stack.push_back(identity_matrix);
    pro_mode = mode;
    if(pro_mode==MGL_PROJECTION)
    {
	 Projection_stack.clear();         
         Projection_stack.push_back(identity_matrix);
        Model_stack.clear();
        Model_stack.push_back(identity_matrix);


    }
#if DEBUG    
    cout<<"the mode style is "<<pro_mode<<endl;
#endif
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
    if(pro_mode == MGL_PROJECTION)
    {
        Projection_stack.push_back(Projection_stack.back());    
    }
    else if (pro_mode==MGL_MODELVIEW)
    {
        Model_stack.push_back(Model_stack.back());    
    }
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
    if(pro_mode == MGL_PROJECTION)
    {
        Projection_stack.pop_back();    
    }
    else if (pro_mode==MGL_MODELVIEW)
    {
        Model_stack.pop_back();    
    }
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
    
    if(pro_mode == MGL_PROJECTION)     
    {
#if DEBUG        
        cout<<"this time is a identiy matrix"<<endl;
#endif       
        Projection_stack.back()=Projection_stack.back()*identity_matrix;
        // Projection_stack.back()=identity_matrix;

    }
    if(pro_mode==MGL_MODELVIEW)
    {
#if DEBUG
        cout<<"load identity for MGL_MODELVIEW"<<endl;
#endif       
        Model_stack.back()=Projection_stack.back()*identity_matrix;    
       //  Model_stack.back()=identity_matrix;    

    }
    
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
    mat4 tmp;
    for(int j=0;j<4;j++)
        for(int i=0;i<4;i++)
            tmp(i,j) = matrix[j*4+i];
    if(pro_mode ==MGL_PROJECTION)
        Projection_stack.back() = tmp;
    else if(pro_mode ==MGL_MODELVIEW)
        Model_stack.back() = tmp;

}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
    mat4 tmp;
    for(int j=0;j<4;j++)
        for(int i=0;i<4;i++)
            tmp(i,j) = matrix[j*4+i];
    if(pro_mode ==MGL_PROJECTION)
        Projection_stack.back() = Projection_stack.back()*tmp;
    else if(pro_mode ==MGL_MODELVIEW)
        Model_stack.back() = Model_stack.back()*tmp;

}


/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
 
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
    //the last column is x,,y,z,1
    mat4 translate_matrix={{1,0,0,0,
                        0,1,0,0,
                        0,0,1,0,
                        x,y,z,1}};
    if(pro_mode ==MGL_PROJECTION)
    {
        Projection_stack.back() = Projection_stack.back()*translate_matrix;     
    }
    if(pro_mode ==MGL_MODELVIEW)
    {
        Model_stack.back() = Model_stack.back()*translate_matrix;     
    }
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
    // angle is 0-360 cos sin use radian
    angle = (angle /180.0f) * PI;

   MGLfloat nor = x*x+y*y+z*z;
    x=x/sqrt(nor);
    y=y/sqrt(nor);
    z=z/sqrt(nor);
   
    mat4 rotate_matrix={{cos(angle)+x*x*(1-cos(angle)),y*x*(1-cos(angle))+z*sin(angle),z*x*(1-cos(angle))-y*sin(angle),0,
                        x*y*(1-cos(angle))-z*sin(angle),cos(angle)+y*y*(1-cos(angle)),z*y*(1-cos(angle))+x*sin(angle),0,
                        x*z*(1-cos(angle))+y*sin(angle),y*z*(1-cos(angle))-x*sin(angle),cos(angle)+z*z*(1-cos(angle)),0,
                        0,0,0,1}};
    
    if(pro_mode ==MGL_PROJECTION)
    {
        Projection_stack.back() = Projection_stack.back()*rotate_matrix;     
    }
    if(pro_mode ==MGL_MODELVIEW)
    {
        Model_stack.back() = Model_stack.back()*rotate_matrix;     
    }
    
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
    mat4 scale_matrix={{x,0,0,0,
                        0,y,0,0,
                        0,0,z,0,
                        0,0,0,1}};
    if(pro_mode ==MGL_PROJECTION)
    {
        Projection_stack.back() = Projection_stack.back()*scale_matrix;     
    }
    if(pro_mode ==MGL_MODELVIEW)
    {
        Model_stack.back() = Model_stack.back()*scale_matrix;     
    }
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{    
    mat4 tmp ={{2.0f*near/(right-left),0.0f,0.0f,0.0f,
                0.0f,2.0f*near/(top-bottom),0.0f,0.0f,
                (right+left)/(right-left),(top+bottom)/(top-bottom),-1*(far+near)/(far-near),-1.0f,
                0.0f,0.0f,(-2.0f*far*near)/(far-near),0.0f}};
    if(pro_mode==MGL_PROJECTION)
    {
#if DEBUG
        printf("make sure this is a projection mode");
#endif       
       // matrix_4 = matrix_4 * tmp;
       Projection_stack.back() = Projection_stack.back()*tmp; 

   }
   if(pro_mode==MGL_MODELVIEW)
   {
#if DEBUG
        cout<<"this time is model_view"<<endl;
#endif       
        Model_stack.back()= Model_stack.back()*tmp;     
   }

}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
    
    mat4 tmp ={{2.0f/(right-left),0.0f,0.0f,0.0f,
                0.0f,2.0f/(top-bottom),0.0f,0.0f,
                0.0f,0.0f,(-2.0f/(far-near)),0.0f,
                -1*(right+left)/(right-left),-(top+bottom)/(top-bottom),-(far+near)/(far-near),1.0f}};
    if(pro_mode==MGL_PROJECTION)
    {
#if DEBUG
        printf("make sure this is a projection mode");
#endif       
       // matrix_4 = matrix_4 * tmp;

    //    Projection_stack.clear();   
             
       Projection_stack.push_back(identity_matrix);
       Projection_stack.back() = Projection_stack.back()*tmp; 

   }
   if(pro_mode==MGL_MODELVIEW)
   {
#if DEBUG
        cout<<"this time is model_view"<<endl;
#endif        
       
      //  Model_stack.clear(); 
        Model_stack.push_back(identity_matrix);
        Model_stack.back()= Model_stack.back()*tmp;     
   }

}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
  //  cout<<"color is set"<<endl;
    RGB[0] = red;
    RGB[1] = green;
    RGB[2] = blue;

}
