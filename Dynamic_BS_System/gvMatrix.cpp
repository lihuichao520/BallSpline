// Project:	VolGrafx
// Group:	VolModel
// File:	coord.cpp
// Date:	05/03/1998
//
// Description:
// routine for coordinates computation
// Designed by Dr. Wu Zhongke
// Include

#include <cmath>

#include "stdafx.h"
#include "gvVector.h"

// Define

#ifndef NULL
#define NULL 0
#endif


gvFLOAT determinate(gvPoint v1,gvPoint v2,gvPoint v3)
{
   return(v1.x*v2.y*v3.z + v1.y*v2.z*v3.x + v1.z*v2.x*v3.y - v1.x*v2.z*v3.y - v1.y*v2.x*v3.z - v1.z*v2.y*v3.x ) ;
}

gvFLOAT det2x2(gvFLOAT a, gvFLOAT b, gvFLOAT c, gvFLOAT d);

/*
Matrix Inversion
by Richard Carling
from "Graphics Gems", Academic Press, 1990
Modified by Zhongke 
*/
gvFLOAT det3x3(gvFLOAT a1,gvFLOAT a2,gvFLOAT a3,gvFLOAT b1,gvFLOAT b2,gvFLOAT b3,gvFLOAT c1,gvFLOAT c2,gvFLOAT c3);

/* 
 *   adjoint( original_matrix, inverse_matrix )
 * 
 *     calculate the adjoint of a 4x4 matrix
 *
 *      Let  a   denote the minor determinant of matrix A obtained by
 *            ij
 *
 *      deleting the ith row and jth column from A.
 *
 *                    i+j
 *     Let  b   = (-1)    a
 *           ij            ji
 *
 *    The matrix B = (b  ) is the adjoint of A
 *                     ij
 */

void adjoint4(gvFLOAT *in, gvFLOAT *out)
{
  gvFLOAT a1, a2, a3, a4, b1, b2, b3, b4;
  gvFLOAT c1, c2, c3, c4, d1, d2, d3, d4;

  /* assign to individual variable names to aid  */
  /* selecting correct values  */
  a1 = in[0]; b1 = in[1]; 
  c1 = in[2]; d1 = in[3];
  a2 = in[4]; b2 = in[5]; 
  c2 = in[6]; d2 = in[7];
  a3 = in[8]; b3 = in[9];
  c3 = in[10]; d3 = in[11];
  a4 = in[12]; b4 = in[13]; 
  c4 = in[14]; d4 = in[15];
  /* row column labeling reversed since we transpose rows & columns */
  out[0]  =   det3x3( b2, b3, b4, c2, c3, c4, d2, d3, d4);
  out[4]  = - det3x3( a2, a3, a4, c2, c3, c4, d2, d3, d4);
  out[8]  =   det3x3( a2, a3, a4, b2, b3, b4, d2, d3, d4);
  out[12]  = - det3x3( a2, a3, a4, b2, b3, b4, c2, c3, c4);

  out[1]  = - det3x3( b1, b3, b4, c1, c3, c4, d1, d3, d4);
  out[5]  =   det3x3( a1, a3, a4, c1, c3, c4, d1, d3, d4);
  out[9]  = - det3x3( a1, a3, a4, b1, b3, b4, d1, d3, d4);
  out[13]  =   det3x3( a1, a3, a4, b1, b3, b4, c1, c3, c4);
      
  out[2]  =   det3x3( b1, b2, b4, c1, c2, c4, d1, d2, d4);
  out[6]  = - det3x3( a1, a2, a4, c1, c2, c4, d1, d2, d4);
  out[10]  =   det3x3( a1, a2, a4, b1, b2, b4, d1, d2, d4);
  out[14]  = - det3x3( a1, a2, a4, b1, b2, b4, c1, c2, c4);
      
  out[3]  = - det3x3( b1, b2, b3, c1, c2, c3, d1, d2, d3);
  out[7]  =   det3x3( a1, a2, a3, c1, c2, c3, d1, d2, d3);
  out[11]  = - det3x3( a1, a2, a3, b1, b2, b3, d1, d2, d3);
  out[15]  =   det3x3( a1, a2, a3, b1, b2, b3, c1, c2, c3);
}
//Matrix3
void adjoint3(gvFLOAT *in, gvFLOAT *out)
{
  gvFLOAT a1, a2, a3, b1, b2, b3;
  gvFLOAT c1, c2, c3;

  /* assign to individual variable names to aid  */
  /* selecting correct values  */
  a1 = in[0]; b1 = in[1]; c1 = in[2]; 
  a2 = in[3]; b2 = in[4]; c2 = in[5]; 
  a3 = in[6]; b3 = in[7]; c3 = in[8]; 
  /* row column labeling reversed since we transpose rows & columns */
  out[0]  =   det2x2( b2, b3, c2, c3);
  out[3]  = - det2x2( a2, a3, c2, c3);
  out[6]  =   det2x2( a2, a3, b2, b3);

  out[1]  = - det2x2( b1, b3, c1, c3);
  out[4]  =   det2x2( a1, a3, c1, c3);
  out[7]  = - det2x2( a1, a3, b1, b3);
      
  out[2]  =   det2x2( b1, b2, c1, c2);
  out[5]  = - det2x2( a1, a2, c1, c2);
  out[8]  =   det2x2( a1, a2, b1, b2);
}


/*
 * gvFLOAT = det2x2( gvFLOAT a, gvFLOAT b, gvFLOAT c, gvFLOAT d )
 * 
 * calculate the determinant of a 2x2 matrix.
 */

gvFLOAT det2x2(gvFLOAT a, gvFLOAT b, gvFLOAT c, gvFLOAT d)
{
  gvFLOAT ans;
  ans = a * d - b * c;
  return ans;
}
/*
 * gvFLOAT = det3x3(  a1, a2, a3, b1, b2, b3, c1, c2, c3 )
 * 
 * calculate the determinant of a 3x3 matrix
 * in the form
 *
 *     | a1,  b1,  c1 |
 *     | a2,  b2,  c2 |
 *     | a3,  b3,  c3 |
 */

gvFLOAT det3x3(gvFLOAT a1,gvFLOAT a2,gvFLOAT a3,gvFLOAT b1,gvFLOAT b2,gvFLOAT b3,gvFLOAT c1,gvFLOAT c2,gvFLOAT c3)
{
  gvFLOAT ans;
  ans = a1 * det2x2( b2, b3, c2, c3 )
      - b1 * det2x2( a2, a3, c2, c3 )
      + c1 * det2x2( a2, a3, b2, b3 );
  return ans;
}

/*
 * gvFLOAT = det4x4( matrix )
 * 
 * calculate the determinant of a 4x4 matrix.
 */
gvFLOAT det4x4(gvFLOAT * m )
{
  gvFLOAT ans;
  gvFLOAT a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3,d4;
  /* assign to individual variable names to aid selecting */
  /*  correct elements */
  a1 = m[0]; b1 = m[1]; c1 = m[2]; d1 = m[3];
  a2 = m[4]; b2 = m[5]; c2 = m[6]; d2 = m[7];
  a3 = m[8]; b3 = m[9]; c3 = m[10]; d3 = m[11];
  a4 = m[12]; b4 = m[13];c4 = m[14]; d4 = m[15];
  ans = a1 * det3x3( b2, b3, b4, c2, c3, c4, d2, d3, d4)
      - b1 * det3x3( a2, a3, a4, c2, c3, c4, d2, d3, d4)
      + c1 * det3x3( a2, a3, a4, b2, b3, b4, d2, d3, d4)
      - d1 * det3x3( a2, a3, a4, b2, b3, b4, c2, c3, c4);
  return ans;
}

//Matrix
void multiplyMatrixmn(gvINT  m, gvINT  n, gvINT  r, gvFLOAT *a, gvFLOAT *b, gvFLOAT *c)
{
	gvINT  i, j, k;

	for(j=0; j<m; j++)
	{
	   for(i=0; i<r; i++)
	   {
		   c[j*r+i] =0.0;
		   for(k=0; k<n; k++)
				c[j*r+i] += a[k*m+j]*b[k*r+i];
	   }
	}
}

void multiplyMatrix(gvINT n, gvFLOAT *a, gvFLOAT *b, gvFLOAT *c)
{
  gvINT  i,j,k;
  
  for(j=0; j<n; j++)
  {
    for(i=0; i<n; i++)
    {
      c[j*n+i] = 0.0;
      for(k=0; k<n; k++)
        c[j*n+i] += a[j*n+k]*b[k*n+i];
    }
  }
}

void multiplyMatrixgeneral(gvINT  m, gvINT  n, gvINT  r, gvFLOAT *a, gvFLOAT *b, gvFLOAT *c)
{
	gvINT  i, j, k;

	for(j=0; j<m; j++)
	{
	   for(i=0; i<r; i++)
	   {
		   c[j*r+i] =0.0;
		   for(k=0; k<n; k++)
				c[j*r+i] += a[j*m+k]*b[k*r+i];
	   }
	}
}

void TransposeMatrix(gvINT  cols, gvINT  rows, gvFLOAT *inM, gvFLOAT *outM)
{
    int tempI, tempJ;
    for(tempI=0; tempI < rows; tempI++)
	for(tempJ=0; tempJ < cols; tempJ++)
	    outM[tempI*cols+tempJ] = inM[tempJ*rows+tempI];    
}

void MultMatrix(gvINT  firstrows, gvINT  cols, gvINT  secondcols, gvFLOAT *firstM, gvFLOAT *secondM, gvFLOAT *outM)
{
    int i,j,k;
    gvFLOAT sum;

    for(i=0; i < secondcols; i++)
	for(j=0; j < firstrows; j++)
	{
	    sum = 0.0;
	    for(k=0; k < cols; k++)
		  sum += firstM[j*cols+k] * secondM[k*secondcols+i];
	    outM[j*secondcols+i] = sum;
	}
}
/* 
 *   inverse( original_matrix, inverse_matrix )
 * 
 *    calculate the inverse of a 4x4 matrix
 *
 *     -1     
 *     A  = ___1__ adjoint A
 *           det A
 */

gvINT  inverse4(gvFLOAT *in,gvFLOAT *out) 
{
  int i, j;
  gvFLOAT det;

  /* calculate the adjoint matrix */
  adjoint4( in, out );
  /*  calculate the 4x4 determinant
   *  if the determinant is zero, 
   *  then the inverse matrix is not unique.
   */
  det = det4x4( in );
  if(fabs( det ) < SMALL_NUMBER)
	return 0;	
    /* scale the adjoint matrix to get the inverse */
  for (i=0; i<4; i++)
    for(j=0; j<4; j++)
      out[i*4+j] = out[i*4+j] / det;
  return 1;
}

gvINT  inverse3(gvFLOAT *in,gvFLOAT *out) 
{
  int i, j;
  gvFLOAT det;

  /* calculate the adjoint matrix */
  adjoint3( in, out );
  /*  calculate the 4x4 determinant
   *  if the determinant is zero, 
   *  then the inverse matrix is not unique.
   */
  det = det3x3(in[0],in[1], in[2],
			   in[3],in[4], in[5],
			   in[6],in[7], in[8]);
  if(fabs( det ) < SMALL_NUMBER)
	return 0;	
    /* scale the adjoint matrix to get the inverse */
  for (i=0; i<3; i++)
    for(j=0; j<3; j++)
      out[i*3+j] = out[i*3+j] / det;
  return 1;
}

gvINT  inversepseudo(gvFLOAT *in, gvFLOAT * out)
{
  gvFLOAT m34[12], m3[9], invm3[9];

  TransposeMatrix(4,3,in,m34); 
  MultMatrix(3,4,3, m34, in, m3);
  if(inverse3(m3, invm3))
  {
	MultMatrix(3,3,4, invm3, m34, out);
	return 1;
  }
  return 0;
}

gvINT  inversepseudo5(gvFLOAT *in, gvFLOAT * out)
{
  gvFLOAT m54[20], m4[16], invm4[16];

  TransposeMatrix(4,5,in,m54); 
  MultMatrix(4,5,4, in, m54, m4);
  if(inverse3(m4, invm4))
  {
	MultMatrix(5,4,4, m54, invm4, out);
	return 1;
  }
  return 0;
}

gvINT  inversepseudo6(gvFLOAT *in, gvFLOAT * out)
{
  gvFLOAT m64[24], m4[16], invm4[16];

  TransposeMatrix(4,6,in,m64); 
  MultMatrix(4,6,4, in, m64, m4);
  if(inverse4(m4, invm4))
  {
	MultMatrix(6,4,4, m64, invm4, out);
	return 1;
  }
  return 0;
}

void multiplyMatrix4(Matrix4 *a, Matrix4 *b, Matrix4 *c)
{
	gvINT  i, j, k;

	for(j=0; j<4; j++)
	{
	   for(i=0; i<4; i++)
	   {
		   c->element[j][i] =0.0;
		   for(k=0; k<4; k++)
				c->element[j][i] += a->element[j][k]*b->element[k][i];
	   }
	}
}

void multiplyMatrix3(Matrix3 *a, Matrix3 *b, Matrix3 *c)
{
	gvINT  i, j, k;

	for(j=0; j<3; j++)
	{
	   for(i=0; i<3; i++)
	   {
		   c->element[j][i] =0.0;
		   for(k=0; k<3; k++)
				c->element[j][i] += a->element[j][k]*b->element[k][i];
	   }
	}
}

