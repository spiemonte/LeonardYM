#include "LinearAlgebraDriver.h"

void dot(double& result_re, double& result_im, const double* ptr1, const double* ptr2, int start, int end) {
	result_re = 0.;
	result_im = 0.;
	for (int index = start; index < end; index += 8) {
		vector4double v1r = {ptr1[index+0],ptr1[index+2],ptr1[index+4],ptr1[index+6]};
		vector4double v1i = {ptr1[index+1],ptr1[index+3],ptr1[index+5],ptr1[index+7]};
		vector4double v2r = {ptr2[index+0],ptr2[index+2],ptr2[index+4],ptr2[index+6]};
		vector4double v2i = {ptr2[index+1],ptr2[index+3],ptr2[index+5],ptr2[index+7]};

		double __attribute__((aligned (32))) tmp[4];

		vec_st(vec_add(vec_mul(v1r,v2r),vec_mul(v1i,v2i)),0L,tmp);
		result_re += tmp[0] + tmp[1] + tmp[2] + tmp[3];
		vec_st(vec_sub(vec_mul(v1r,v2i),vec_mul(v1i,v2r)),0L,tmp);
		result_im += tmp[0] + tmp[1] + tmp[2] + tmp[3];
	}
}
