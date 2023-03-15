This is the two-dimensional MNOMP algorihtm for the following paper
L. Han, X. Liu, N. Zhang, S. Wu, J. Zhu and Z. Xu,  Two-dimensional multi-snapshot Newtonalized orthogonal matching pursuit for DOA estimation, Digital Signal Processing, vol. 121, Mar. 2022.
The code is written by Jiang Zhu, Menghuai Xu and Lin Han. If you have any problems, please feel free to contact me with jiangzhu16@zju.edu.cn
The following files are:
1. demo_MNOMP2D.m: the main function to test 2D-NOMP algorihtm.
2. NOMP2D.m: the 2D-NOMP function
3. DectectNew_2D: to detect a new target from the 2D line spectrum
4. LeastSquare_2D: to adjuste the complex amplitudes of targets by least
    square method
5. RefineOne_2D.m: to refine the angular frequency of a particular target with
    NOMP method
6. RefineAll_2D.m: to refine the angular frequency of all the targets one by
    one.




