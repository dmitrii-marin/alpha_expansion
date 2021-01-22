
typedef float Dtype;

void run_expansion(const Dtype *images, int, int, int, int,
                   const Dtype *init_sol, int, int, int, int, const Dtype *ROI,
                   int, int, int, int, const Dtype *seeds, int, int, int, int,
                   const Dtype *unary, int N, int C, int H, int W, int max_iter,
                   float potts_weight, Dtype *gc_segmentation, int, int, int,
                   int, Dtype *unary_energy, int, Dtype *smooth_energy, int);
