//****************************************************************************//
//                         P.cpp                                              //
//                         description: the implementation of P.cpp           //
//                         function   : Control Points and Control Radius     //
//                                      Vector P(t)                           //
//                         author     : lihuichao                             //
//                         Date       : 2016/11/24                            //
//****************************************************************************//
#include "stdafx.h"
#include "P.h"


P::P()
{
}


P::~P()
{
}

void P::InitialP(int size){

	Size = size;
	Pd = new double[size];
	Pt = new double[size];
	Ptt = new double[size];

	ResetZeroP();
	ResetZeroPt();
	ResetZeroPtt();
}

void P::Destroy(){

	delete[] Pd;
	delete[] Pt;
	delete[] Ptt;
}

void P::ResetZeroP(){

	for (int i = 0; i < Size; i++){
		Pd[i] = 0;
	}
}

void P::ResetZeroPt(){
	for (int i = 0; i < Size; i++){
		Pt[i] = 0;
	}
}

void P::ResetZeroPtt(){
	for (int i = 0; i < Size; i++){
		Ptt[i] = 0;
	}
}