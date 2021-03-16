#include "alpha_expansion.hpp"

#include "gco-v3.0/GCoptimization.h"

#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>
#include <iostream>

#include <float.h>
#include <stdio.h>

#include <omp.h>

static const int UNLABELED = 255;
static const int dh[] = { 0, 1, 1, -1 };
static const int dw[] = { 1, 0, 1,  1 };

void run_expansion_single(const Dtype *image, int K, const Dtype *init_sol,
                          const Dtype *ROI, const Dtype *seeds,
                          const Dtype *unary, int C, int H, int W, int max_iter,
                          float potts_weight, Dtype *gc_segmentation,
                          Dtype *unary_energy, Dtype *smooth_energy) {

  const double SCALE = 10000. / potts_weight;
  const double HARD = 10 * potts_weight;

  std::unique_ptr<GCoptimizationGeneralGraph> gc(
      new GCoptimizationGeneralGraph(W * H, C));

  double time = omp_get_wtime();
  // printf("%.2lf: Entered\n", omp_get_wtime() - time);

  if (init_sol)
    for (int h = 0; h < H; h++)
      for (int w = 0; w < W; w++)
        if (ROI[h * W + w]) {
          gc->setLabel(h * W + w, init_sol[h * W + w]);
        }

  // for(int h=0; h<H; h++) {
  //   for(int w=0; w<W; w++) {
  //     if (ROI[h*W + w]==false) {
  //       gc_segmentation[h*W + w] = 0;
  //       continue;
  //     }
  //     Dtype minsofar = INFINITY;
  //     Dtype label = -1;
  //     Dtype seed_i = seeds[h*W + w];
  //     if (seed_i != UNLABELED) {
  //       label = seed_i;
  //     } else {
  //       for(int c=0; c<C; c++) {
  //         Dtype val = unary[c*H*W + h*W + w];
  //         if (val != 0 && !std::isnormal(val)) {
  //           throw std::runtime_error("unaries contain invalid values");
  //         }
  //         if (val < minsofar){
  //           minsofar = unary[c*H*W + h*W + w];
  //           label = c;
  //         }
  //       }
  //     }
  //     gc->setLabel(h*W + w, label);
  //   }
  // }

  // first set up data costs individually
  int count = 0;
  for (int c = 0; c < C; c++) {
    for (int h = 0; h < H; h++) {
      for (int w = 0; w < W; w++) {
        if (ROI[h * W + w] == false)
          continue;
        Dtype seed_i = seeds[h * W + w];
        Dtype cost;
        if (seed_i == UNLABELED)
          cost = unary[c * H * W + h * W + w];
        else {
          if (int(seed_i) >= C || int(seed_i) < 0) {
            fprintf(stderr, "%d\n", int(seed_i));
            throw std::runtime_error("bad value");
          }
          cost = unary[int(seed_i) * H * W + h * W + w];
          if (seed_i != c)
            cost += HARD;
        }
        gc->setDataCost(h * W + w, c, SCALE * cost);
        ++count;
      }
    }
  }
  count /= C;

  // printf("%.2lf: Unaries are set\n", omp_get_wtime() - time);

  // next set up the array for smooth costs
  std::vector<GCoptimizationGeneralGraph::EnergyTermType> smooth(C * C);
  for (int l1 = 0; l1 < C; l1++)
    for (int l2 = 0; l2 < C; l2++)
      smooth[l1 + l2 * C] = (l1 != l2) ? 1 : 0;
  gc->setSmoothCost(smooth.data());


  // expected squared distance
  double diff_total = 0;
  int diff_count = 0;
  for (int c = 0; c < K; ++c) {
    for (int i = 0; i < 4; i++) {
      for (int h = 0; h < H; h++) {
	for (int w = 0; w < W; w++) {
	  if (ROI[h * W + w] == false)
	    continue;
	  int h2 = h + dh[i];
	  int w2 = w + dw[i];
	  if (w2 < 0 || w2 >= W || h2 < 0 || h2 >= H)
	    continue;
	  if (ROI[h2 * W + w2] == false)
	    continue;
	  if (K == 3) {
		  Dtype diff =
		      image[c * H * W + h * W + w] - image[c * H * W + h2 * W + w2];
		  diff_total += diff * diff;
	  }
	  diff_count++;
	}
      }
    }
  }
  diff_count /= K;
  double sigmasquared = diff_total / diff_count;

  if (sigmasquared == 0)
    sigmasquared = 1;

  // printf("%.2lf: pairwise sigma is done\n", omp_get_wtime() - time);
  // printf("sigma2 is %f (%d)\n", sigmasquared * 2, diff_count);

  // now set up a grid neighborhood system
  for (int h = 0; h < H; h++)
    for (int w = 0; w < W; w++) {
      if (ROI[h * W + w] == false)
        continue;
      for (int i = 0; i < 4; i++) {
        int h2 = h + dh[i];
        int w2 = w + dw[i];
        if (w2 < 0 || w2 >= W || h2 < 0 || h2 >= H)
          continue;
        if (ROI[h2 * W + w2] == false)
          continue;
        Dtype cost;
        if (K == 3) {
		Dtype diff_b =
		    image[0 * H * W + h * W + w] - image[0 * H * W + h2 * W + w2];
		Dtype diff_g =
		    image[1 * H * W + h * W + w] - image[1 * H * W + h2 * W + w2];
		Dtype diff_r =
		    image[2 * H * W + h * W + w] - image[2 * H * W + h2 * W + w2];
		Dtype diff = diff_b * diff_b + diff_g * diff_g + diff_r * diff_r;
		cost = potts_weight * exp(-diff / 2 / sigmasquared);
		if (i > 1)
		  cost /= M_SQRT2;
		cost = cost * count / diff_count;
        } else {
		Dtype diff = std::max(image[h * W + w], image[h2 * W + w2]);
		cost = potts_weight * std::max(Dtype(0), 1 - diff) * count / diff_count;
        }
	gc->setNeighbors(h * W + w, h2 * W + w2, SCALE * cost);
      }
    }

  // printf("%.2lf: pairwise set\n", omp_get_wtime() - time);

  // LOG(INFO) << "sigmasquared is " << sigmasquared;

  // printf("\nBefore optimization energy is %d\n",gc->compute_energy());
  // printf("max_iter: %d\n",max_iter);
  // try {
  gc->expansion(max_iter); // run expansion for max_iter. For swap use
                           // gc->swap(num_iterations);
  // } catch (const GCException& ex) {
  //   LOG(FATAL) << "GCException is thrown" << std::endl << ex.message;
  // }
  // LOG(INFO) << "After optimization energy is " << gc->compute_energy() /
  // double(SCALE);

  // printf("%.2lf: alpha done\n", omp_get_wtime() - time);

  for (int h = 0; h < H; h++) {
    for (int w = 0; w < W; w++) {
      if (ROI[h * W + w] == false)
        continue;

      Dtype label = gc_segmentation[h * W + w] = gc->whatLabel(h * W + w);
#if 0
      Dtype seed_i = seeds[h * W + w];
      if (seed_i != UNLABELED && seed_i != label) {
        std::cerr << "Found label " << label << " at (" << h << ", " << w << ")"
          << " while expecting " << seed_i
          << " [init_sol: " << init_sol[h * W + w] << "]" << std::endl;
        throw std::runtime_error("Failed to respect the hard constraints");
      }
#endif
    }
  }
  if (unary_energy)
    *unary_energy = Dtype(gc->giveDataEnergy()) / SCALE / count;
  if (smooth_energy)
    *smooth_energy = Dtype(gc->giveSmoothEnergy()) / SCALE / count;
  // printf("%lf %d %lf\n", *unary_energy, count, *smooth_energy);
}

void run_expansion_batch(const Dtype *image, int K, const Dtype *init_sol,
                         const Dtype *ROI, const Dtype *seeds,
                         const Dtype *unary, int N, int C, int H, int W,
                         int max_iter, float potts_weight,
                         Dtype *gc_segmentation,
                         Dtype *unary_energy,
                         Dtype *smooth_energy)
{
#pragma omp parallel for
  for (int i = 0; i < N; ++i)
    try {
      run_expansion_single(image + i * K * H * W, K, init_sol ? init_sol + i * H * W : NULL,
                           ROI + i * H * W, seeds + i * H * W,
                           unary + i * C * H * W, C, H, W, max_iter, potts_weight,
                           gc_segmentation + i * H * W, unary_energy + i,
                           smooth_energy + i);
    } catch (GCException& ex) {
      ex.Report();
      throw;
    }
}

void run_expansion(const Dtype *images, int, int K, int, int,
                   const Dtype *init_sol, int, int, int, int, const Dtype *ROI,
                   int, int, int, int, const Dtype *seeds, int, int, int, int,
                   const Dtype *unary, int N, int C, int H, int W, int max_iter,
                   float potts_weight, Dtype *gc_segmentation, int, int, int,
                   int, Dtype *unary_energy, int, Dtype *smooth_energy, int)
{
  run_expansion_batch(images, K, init_sol, ROI, seeds, unary, N, C, H, W, max_iter,
                      potts_weight, gc_segmentation, unary_energy, smooth_energy);
}
