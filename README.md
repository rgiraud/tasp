This code is free to use for any non-commercial purposes.
It contains an implementation of the TASP superpixel method proposed in:
[1] - Texture-Aware Superpixel Segmentation
      Rémi Giraud, Vinh-Thong Ta, Nicolas Papadakis, Yannick Berthoumieu
      IEEE International Conference on Image Processing (ICIP 2019)

Note that the core of the implementation is based on the provided code associated to the following paper:
[2] - Zhengqin Li, Jiansheng Chen
      Superpixel Segmentation using Linear Spectral Clustering
      International Conference on Computer Vision and Pattern Recognition (CVPR), 2015

If you use this code, please cite both [1] and [2].

(C) Rémi Giraud, 2019
rgiraud@u-bordeaux.fr, https://remi-giraud.enseirb-matmeca.fr
Bordeaux-INP, IMS Laboratory


%%%%%%% MATLAB - C++/MEX %%%%%%%%%

run main.m



%%%%%%% C++ %%%%%%%%%

////// COMPILATION /////

- Need CImg library :  http://cimg.eu/

- To compile : make


////// EXECUTION //////

./TASP -i img_name [-outm output_map_name] [-outb output_border_name] [-k superpixel_nbr] [-m compactness] [-Kpm knn_by_PM] [-patch_w patch_size] 

Example with contour map:  (make test)
./TASP -i test_img.jpg -k 300 -m 0.1 -outm test_img_map.png -outb test_img_border.png

On an image list:  (make test_list)
./scripts/test_list.sh ./data/list_file.txt ./data/ 450 0.1



