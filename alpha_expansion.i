%module alphaexpansion


%{
#define SWIG_FILE_WITH_INIT
#include "alpha_expansion.hpp"
%}

%include "../bilateralfilter/numpy.i"

%init %{
    import_array();
%}

%apply (float* IN_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4) { (const Dtype* images, int, int, int, int) }
%apply (float* IN_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4) { (const Dtype* init_sol, int, int, int, int) }
%apply (float* IN_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4) { (const Dtype* ROI, int, int, int, int) }
%apply (float* IN_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4) { (const Dtype* seeds, int, int, int, int) }
%apply (float* IN_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4) { (const Dtype* unary, int N, int C, int H, int W) }
%apply (float* INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4) { (Dtype* gc_segmentation, int, int, int, int) }
%apply (float* INPLACE_ARRAY1, int DIM1) { (Dtype* unary_energy, int), (Dtype* smooth_energy, int) }

%include "alpha_expansion.hpp"
