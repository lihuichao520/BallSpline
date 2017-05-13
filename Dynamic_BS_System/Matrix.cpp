//****************************************************************************//
//                         Matrix.cpp                                         //
//                         description: the implementation of Matrix.cpp      //
//                         function   : Matrix(mass,damp,stiffness)           //
//                         author     : lihuichao                             //
//                         Date       : 2016/11/22                            //
//****************************************************************************//
#include "stdafx.h"
#include "Matrix.h"


Matrix::Matrix()
{
}


Matrix::~Matrix()
{
}


void Matrix::InitialMatrix(int row,int column){
	
	M = new double *[row];

	for (int i = 0; i < row; i++){

		M[i] = new double[column];
	}

	Row = row;
	Column = column;
	ResetZero();
}
void Matrix::DestroyMatrix(){
	for (int i = 0; i < Row; i++){

		delete[] M[i];
	}

	delete M;
}
void Matrix::ResetZero(){

	for (int i = 0; i < Row; i++){

		for (int j = 0; j < Column; j++){

			M[i][j] = 0;
		}
	}
}