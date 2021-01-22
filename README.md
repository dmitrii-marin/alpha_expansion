# Alpha Expansion

This is a multi-thread wrapper over [GCO v3](https://vision.cs.uwaterloo.ca/code/). 

## Compile and install

```sh
# swig -python -c++ alpha_expansion.i
# python setup.py install
```

## Use

```python

#  alphaexpansion.run_expansion runs alpha expansion algorithm on images `images`
#
#  Parameters:
#       `images` - a numpy array (float32) of size (N, C, H, W) where N is the number of images, 
#                  C is the number of color channels in the image (C=3), H and W are the 
#                  dimensions of each image
#       `x0`     - a numpy array (float32) of size (N, 1, H, W) with each value in {0, 1, .., K-1}   
#                  where K is the number of clusters
#       `ROI`    - a numpy array (float32) of size (N, 1, H, W) with each value either 0 or 1
#                  where 0 indicated pixels that should not participate in the optimization. 
#       `seeds`  - a numpy array (float32) of size (N, 1, H, W) with each value in {0, 1, .., K-1, 255}
#                  where values below K indicate pixels with fixed labels
#       `unary`  - a numpy array (float32) of size (N, K, H, W) with unary potentials for each pixel
#       `max_iter` - number of alpha-expansion iterations (one iterations goes over all labels once)
#       `potts_weight` - is the weight of Potts term in alpha exapnsion
#
# Returns:
#       `out`            - a numpy array (float32) of size (N, 1, H, W) with the result labels
#       `unary_energy`   - a numpy array (float32) of size (N) with unary energy value at the 
#                          end of the iterations
#       `smooth_energy`  - a numpy array (float32) of size (N) with pair-wise energy value at the 
#                          end of the iterations
#
def alphaexpansion.run_expansion(
            images, x0, ROI, seeds, unary,
            max_iter, potts_weight, out, unary_energy, smooth_energy)
```

## Citation

If you find this code useful in your research, consider citing
```
@InProceedings{ADM:cvpr19,
  author = {Dmitrii Marin and Meng Tang and Ismail Ben Ayed and Yuri Boykov},
  title = {Beyond Gradient Descent for Regularized Segmentation Losses},
  booktitle = {IEEE conference on Computer Vision and Pattern Recognition (CVPR)},
  volume = {},
  pages = {},
  month = {June},
  year = {2019},
  address = {Long Beach, California}
}
```
Also read the underlying code's [README](gco-v3.0/GCO_README.TXT)
