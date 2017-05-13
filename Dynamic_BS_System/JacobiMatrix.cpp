//****************************************************************************//
//                         JacobiMatrix.cpp                                   //
//                         description: the implementation of JacobiMatrix.cpp//
//                         function   : Jacobi Matrix(time)                   //
//                         author     : lihuichao                             //
//                         Date       : 2016/11/22                            //
//****************************************************************************//

#include "stdafx.h"
#include "JacobiMatrix.h"



JacobiMatrix::JacobiMatrix()
{
}

JacobiMatrix::~JacobiMatrix()
{
}

void JacobiMatrix::Initial(int row, int column){

	Row = row;
	Column = column;

	//lhc-- ����ռ䲢��ʼ��Ϊ0
	J1.zeros(row, column);
	Ju1.zeros(row, column);
	Juu1.zeros(row, column);

}

//lhc-- ���ݲ���ֵu�������Ĵ���k�����ſ˱Ⱦ���,������Nik��u����i=0��...��n��
void JacobiMatrix::setJMatrix(double u, int n, int k, DisplayBBSC bbsc, int leftPoint){
	int i;
	int row;
	int column;
	double Nik = 0.;

	//lhc-- ����
	J1.zeros();
	
	//lhc-- ���ݲ���������B����������
	for (i = leftPoint - bbsc.curve.order + 1, column = 0; i <= leftPoint; i++){

		Nik = CalBasisFunc(u, i, k, bbsc);
		for (row = 0; row < Row && column < Column; row++, column++){

			J1(row,column) = Nik;
		}
	}

}

//lhc-- ����һ�׵����ſ˱Ⱦ���Ju,������Nik'��u����i=0��...��n��
void JacobiMatrix::setJuMatrix(double u, int n, int k, DisplayBBSC bbsc, int leftPoint){

	int i;
	int row;
	int column;
	double dNik = 0.;

	//lhc-- ����
	Ju1.zeros();

	//lhc-- ���ݲ���������B����������
	for (i = leftPoint - bbsc.curve.order + 1, column = 0; i <= leftPoint; i++){

		dNik = FirDerOfBasis(u, i, k, bbsc);
		for (row = 0; row < Row && column < Column; row++, column++){

			Ju1(row, column) = dNik;
		}
	}
}

//lhc-- ���ö��׵����ſ˱Ⱦ���Juu��������Nik'(u)(i=0,...,n)
void JacobiMatrix::setJuuMatrix(double u, int n, int k, DisplayBBSC bbsc, int leftPoint){

	int i;
	int row;
	int column;
	double duuNik = 0.;

	//lhc-- ����
	Juu1.zeros();

	//lhc-- ���ݲ���������B����������
	for (i = leftPoint - bbsc.curve.order + 1, column = 0; i <= leftPoint; i++){

		duuNik = SecDerOfBasis(u, i, k, bbsc);
		for (row = 0; row < Row && column < Column; row++, column++){

			Juu1(row, column) = duuNik;
		}
	}
}

//lhc-- ���ܴ�������(2016/12/20û������)
double JacobiMatrix::CalBasisFunc(double u, int i, int k, DisplayBBSC bbsc){

	double val1 = 0;
	double val2 = 0;
	double Val = 0.0;

	if (k == 0){

		if (u < bbsc.curve.knot[i] || u > bbsc.curve.knot[i+1] ){
			return Val;
		}
		else{

			Val = 1.0;
			return Val;
		}
	}

	if (k > 0){
		if ( u < bbsc.curve.knot[i] || u > bbsc.curve.knot[i + k + 1]){
			return Val;
		}
		else{

			double alpha = 0.0;
			double beta = 0.0;
			double dTemp = 0.0;

			dTemp = bbsc.curve.knot[i + k] - bbsc.curve.knot[i];

			if (dTemp == 0.0){

				alpha = 0.0;
			}
			else{

				alpha = (u - bbsc.curve.knot[i]) / dTemp;
			}

			dTemp = bbsc.curve.knot[i + k + 1] - bbsc.curve.knot[i + 1];

			if (dTemp == 0.0){

				beta = 0.;
			}
			else{
				beta = (bbsc.curve.knot[i + k + 1] - u) / dTemp;
			}

			val1 = alpha * CalBasisFunc(u, i, k - 1, bbsc);
			val2 = beta * CalBasisFunc(u, i + 1, k - 1, bbsc);
			Val = val1 + val2;
		}// end of else
	}// end of if

	return Val;
}

//lhc-- ���ݲ���ֵ�ʹ���k�������i����������һ�׵���
double JacobiMatrix::FirDerOfBasis(double u, int i, int k, DisplayBBSC bbsc){

	double val1 = 0;
	double val2 = 0;
	double Val = 0.0;
	double alpha = 0.0;
	double beta = 0.0;
	double dTemp = 0.0;

	dTemp = bbsc.curve.knot[i + k] - bbsc.curve.knot[i];
	if (dTemp == 0.0){

		alpha = 0.0;
	}
	else{

		alpha = 1 / dTemp;
	}

	dTemp = bbsc.curve.knot[i + k + 1] - bbsc.curve.knot[i + 1];
	if (dTemp == 0.0){

		beta = 0.;
	}
	else{
		beta = 1 / dTemp;
	}

	val1 = alpha * CalBasisFunc(u, i, k - 1, bbsc);
	val2 = beta * CalBasisFunc(u, i + 1, k - 1, bbsc);
	Val = k*(val1 - val2);

	return Val;
}

//lhc-- ���ݲ���ֵ�ʹ���k�������i���������Ķ��׵���
double JacobiMatrix::SecDerOfBasis(double u, int i, int k, DisplayBBSC bbsc){

	double val1 = 0;
	double val2 = 0;
	double Val = 0.0;
	double alpha = 0.0;
	double beta = 0.0;
	double dTemp = 0.0;

	dTemp = bbsc.curve.knot[i + k] - bbsc.curve.knot[i];
	if (dTemp == 0.0){

		alpha = 0.0;
	}
	else{

		alpha = 1 / dTemp;
	}

	dTemp = bbsc.curve.knot[i + k + 1] - bbsc.curve.knot[i + 1];
	if (dTemp == 0.0){

		beta = 0.;
	}
	else{
		beta = 1 / dTemp;
	}

	val1 = alpha * FirDerOfBasis(u, i, k - 1, bbsc);
	val2 = beta * FirDerOfBasis(u, i + 1, k - 1, bbsc);
	Val = k*(val1 - val2);

	return Val;
}



//lhc-- ԭʼ����++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
void JacobiMatrix::InitialJacobiMatrix(int column){

	J = new double *[4];
	for (int i = 0; i < 4; i++){

		J[i] = new double[column];
	}

	Ju = new double *[4];
	for (int i = 0; i < 4; i++){

		Ju[i] = new double[column];
	}

	Juu = new double *[4];
	for (int i = 0; i < 4; i++){

		Juu[i] = new double[column];
	}

	Row = 4;
	Column = column;
	ResetZeroJ();
	ResetZeroJu();
	ResetZeroJuu();
}

//lhc-- ���ݲ���ֵu�������Ĵ���k�����ſ˱Ⱦ���
void JacobiMatrix::SetJacobiMatrix(double u, int n, int k, DisplayBBSC bbsc, int leftPoint){

	int i;
	int row;
	int column;
	double Nik = 0.;


	//lhc-- ����
	ResetZeroJ();
	//lhc-- ���ݲ���������B����������
	for (i = leftPoint - bbsc.curve.order + 1, column = 0; i <= leftPoint; i++){

		Nik = CalBasisFunc(u, i, k, bbsc);
		for (row = 0; row < Row && column < Column; row++, column++){

			J[row][column] = Nik;
		}
	}

}

//lhc-- ���ö�һ�׵����ſ˱Ⱦ���Ju
void JacobiMatrix::SetFirDerJacobiMatrix(double u, int n, int k, DisplayBBSC bbsc, int leftPoint){

	int i;
	int row;
	int column;
	double dNik = 0.;

	//lhc-- ����
	ResetZeroJu();
	//lhc-- ���ݲ���������B����������
	for (i = leftPoint - bbsc.curve.order + 1, column = 0; i <= leftPoint; i++){

		dNik = FirDerOfBasis(u, i, k, bbsc);
		for (row = 0; row < Row && column < Column; row++, column++){

			Ju[row][column] = dNik;
		}
	}
}

//lhc-- ���öԶ��׵����ſ˱Ⱦ���Juu
void JacobiMatrix::SetSecDerJacobiMatrix(double u, int n, int k, DisplayBBSC bbsc, int leftPoint){

	int i;
	int row;
	int column;
	double duuNik = 0.;

	//lhc-- ����
	ResetZeroJuu();
	//lhc-- ���ݲ���������B����������
	for (i = leftPoint - bbsc.curve.order + 1, column = 0; i <= leftPoint; i++){

		duuNik = SecDerOfBasis(u, i, k, bbsc);
		for (row = 0; row < Row && column < Column; row++, column++){

			Juu[row][column] = duuNik;
		}
	}
}

void JacobiMatrix::ResetZeroJ(){

	for (int i = 0; i < Row; i++){

		for (int j = 0; j < Column; j++){

			J[i][j] = 0;
		}
	}
}

void JacobiMatrix::ResetZeroJu(){
	for (int i = 0; i < Row; i++){

		for (int j = 0; j < Column; j++){
		
			Ju[i][j] = 0;
		}
	}
}

void JacobiMatrix::ResetZeroJuu(){
	for (int i = 0; i < Row; i++){

		for (int j = 0; j < Column; j++){

			Juu[i][j] = 0;
		}
	}
}

void JacobiMatrix::DestroyJacobiMatrix(){

	for (int i = 0; i < 4; i++){

		delete[] J[i];
		delete[] Ju[i];
		delete[] Juu[i];
	}

	delete J;
	delete Ju;
	delete Juu;
}
//lhc-- ԭʼ����++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

void JacobiMatrix::Destory(){
}