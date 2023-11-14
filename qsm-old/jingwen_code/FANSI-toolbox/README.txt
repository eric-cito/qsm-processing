The FAst Nonlinear Susceptibility Inversion (FANSI) Toolbox was developed by 
Carlos Milovic, PhD. at the Biomedical Imaging Center at Pontificia Universidad 
Catolica de Chile and the Wellcome Trust Centre for Neuroimaging at University 
College London, in 2017, in collaboration with:
Berkin Bilgic, PhD. and Bo Zhao, PhD. at Martinos Center for Biomedical Imaging, 
Harvard Medical School, MA, USA
Julio Acosta-Cabronero, PhD at Wellcome Trust Centre for Neuroimaging, Institute 
of Neurology, University College London, London, UK, and German Center for 
Neurodegenerative Diseases (DZNE), Magdeburg, Germany
and Cristian Tejos, PhD. at Department of Electrical Engineering, Pontificia 
Universidad Catolica de Chile, Santiago, Chile and the Biomedical Imaging Center 
at Pontificia Universidad Catolica de Chile.

From 2020, this toolbox is being mantained and expanded at the Department of 
Medical Physics and Biomedical Engineering, University College London, UK, by 
Carlos Milovic, PhD. with supervision by Karin Shmueli, PhD.


This README file corresponds to the FANSI-Toolbox v2.0, released on 2020.07.27


Source code for the Fast Nonlinear Susceptibility Inversion (FANSI) method is
provided in the root of this toolbox, along with some example scripts.
Please see the scripts to learn how to use the functions included in this 
toolbox. You may also use the "doc" and "help" commands to see the header of 
each function and learn about all the function parameters.

The code in this toolbox is based on source code released by Berkin Bilgic at 
http://martinos.org/~berkin/software.html

Please cite the following publications if you use any part of this toolbox:
- Bilgic B., Chatnuntawech I., Langkammer C., Setsompop K.; Sparse Methods 
for Quantitative Susceptibility Mapping; Wavelets and Sparsity XVI, 
SPIE 2015
- Milovic C, Bilgic B, Zhao B, Acosta-Cabronero J, Tejos C. Fast Nonlinear 
Susceptibility Inversion with Variational Regularization. Magn Reson Med.
Accepted Dec. 12th 2017. DOI: 10.1002/mrm.27073


If you use the frequency-discretized dipole formulation, please cite:
- Carlos Milovic, Julio Acosta-Cabronero, Jose Miguel Pinto, Hendrik Mattern, 
Marcelo Andia, Sergio Uribe, and Cristian Tejos. A new discrete dipole kernel 
for quantitative susceptibility mapping; Magnetic Resonance Imaging, 2018;51:7-13.


While the focus of this toolbox is to provide state-of-the-art methods to perform 
the dipole inversion step in quantitative susceptibility mapping (QSM), we have 
included code to perform tasks throughout  the QSM pipeline for completeness and 
for academic and research purposes. Source code for other methods that extend the 
scope of the toolbox is included in different folders, organized into categories, 
as described below:

+ AnalyticModels: This folder contains code to create spheres and cylinders, and 
        simulate the magnetization field that they create. For spheres, there are 
        two functions that mimic the MRI acquisition by calculating intra-voxel 
        dephasing effects. Please cite: Marques JPP, Bowtell R. Application of a 
        fourier-based method for rapid calculation of field inhomogeneity due to 
        spatial variation of magnetic susceptibility. Concepts Magn Reson Part B 
        Magn Reson Eng. 2005;25:65-78. 
                  
+ BackgroundFieldRemoval: Files in this folder include a function to fit a 
        polynomial to the data inside a region of interest (ROI), a function to 
        generate the magnetization field given all air-tissue interfaces, and a 
        mask to fit a background model with additional linear gradients to the 
        phase data via least squared minimization.
        These two folders contain the following methods:
        - MSMV: Multiscale Spherical Mean Value. An efficient alternative 
                background field removal technique based on a Laplacian pyramid 
                decomposition of VSHARP. Please cite:
                Milovic C, Langkammer C, Uribe S, Irarrazaval P, Acosta-Cabronero 
                J, Tejos C. Multiscale Spherical Mean Value based background 
                field removal method for Quantitative Susceptibility Mapping. 
                27th International Conference of the ISMRM, Montreal, Canada, 2019.
        - NPDF: ADMM based solver for the Projection onto Dipole Fields method, 
                including a nonlinear method, with spatial weighting and 
                regularization. Please cite: Milovic C, Bilgic B, Zhao B, 
                Langkammer C, Tejos C and Acosta-Cabronero J. Weak-harmonic 
                regularization for quantitative susceptibility mapping (WH-QSM). 
                Magn Reson Med, 2019;81:1399-1411.
                     
+ Inversion: Algorithms that perform the final dipole inversion step in the QSM 
        pipeline, starting with an estimate of the local magnetization fields. 
        This folder contains the following folders and methods:
        - Closed Form: Truncated K-space Division (TKD) method, based on Shmueli 
                et al MRM 2009, Direct Tikhonov reconstruction and the L2-norm 
                direct solver (this is actually a Tikhonov term imposed on the 
                Gradients of the solutions), based on Bilgic SPIE 2015.
        - eNDI: Based on Polak et al NMR Biomed 2020 (doi:10.1002/nbm.4271), 
                here you may find both a Gradient Descent and a Conjugate 
                Gradient method to perform nonregularized (efficient) Nonlinear 
                Dipole Inversions. Please cite: Milovic C, Karsa A, Shmueli K. 
                Efficient Early Stopping algorithm for Quantitative Susceptibility 
                Mapping (QSM). ESMRMB virtual meeting 2020.
        - PI-QSM: An L1-norm data fidelity term based QSM method that rejects 
                phase inconsistencies, noise and errors, preventing streaking 
                artifacts. Please cite: Milovic C, Lambert M, Langkammer C, 
                Bredies K, Tejos C, and Irarrazaval P. QSM streaking suppression 
                with L1 data fidelity terms. 28th International Virtual ISMRM
                Conference Meeting, 2020. P3257.
                Milovic C, Lambert M, Langkammer C, Bredies K, Irarrazaval P, 
                and Tejos C. Robust to Phase Inconsistencies Quantitative 
                Susceptibility Mapping (PI-QSM). MRM. Submitted 2020.
        - WH-QSM: Weak harmonic QSM method to jointly reconstruct susceptibility 
                maps and background (harmonic) field remnants. Please cite:
                Milovic C, Bilgic B, Zhao B, Langkammer C, Tejos C and 
                Acosta-Cabronero J. Weak-harmonic regularization for quantitative 
                susceptibility mapping (WH-QSM); Magn Reson Med, 2019;81:1399-1411.
             
+ Metrics: Compute different global quality/error metrics to compare reconstructed 
        maps with known ground-truths. The functions used to calculate quality 
        metrics are based on the source code provided at the 4th International 
        Workshop on MRI Phase Contrast & Quantitative Susceptibility Mapping, 
        September 26th - 28th 2016, at Medical University of Graz, Austria, for 
        its QSM Challenge. http://www.neuroimaging.at/pages/qsm.php 
        If you use the HFEN or SSIM metrics, please cite:
        Langkammer, C., Schweser, F., Shmueli, K., Kames, C., Li, X., Guo, L., 
        Milovic, C., Kim, J., Wei, H., Bredies, K., Buch, S., Guo, Y., Liu, Z., 
        Meineke, J., Rauscher, A., Marques, J. P. and Bilgic, B. (2017), 
        Quantitative susceptibility mapping: Report from the 2016 reconstruction 
        challenge. Magn. Reson. Med. doi:10.1002/mrm.26830
        The function to calculate the SSIM quality index uses code by Zhou Wang. 
        Please see the header information for more details. This code was modified 
        by Berkin Bilgic for the 2016 QSM Challenge. Please cite:
        Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli. Image quality 
        assessment: From error measurement to structural similarity. IEEE 
        Transactions on Image Processing, 2004;13(4):600-612.
        If you use the code designed for the XSIM metric, please cite: Milovic C, 
        Tejos C, and Irarrazaval P. Structural Similarity Index Metric setup for 
        QSM applications (XSIM). 5th International Workshop on MRI Phase Contrast 
        & Quantitative Susceptibility Mapping, Seoul, Korea, 2019. 
        (A paper similarly titled is in preparation).
        The function to calculate the Mutual Information quality index uses code 
        by R. Moddemeijer, at http://www.cs.rug.nl/~rudy/matlab/. Please cite: 
        Moddemeijer, R. On Estimation of Entropy and Mutual Information of 
        Continuous Distributions, Signal Processing, 1989;16(3):233-246
        If you use any of the other metrics, please cite: Milovic C, Tejos C, 
        Acosta-Cabronero J, Ozbay PS, Schweser F, Marques JP, Irarrazaval P, 
        Bilgic B, Langkammer C. The 2016 QSM Challenge: Lessons learned and 
        considerations for a future challenge design; Magn Reson Med. 
        2020;85:1624-1637 doi: 10.1002/mrm.28185

+ MultiechoFit: Nonlinear methods to estimate the field map and phase offset 
        from multiecho GRE acquisitions. Regularization terms are imposed on the 
        phase offset to reduce noise.

+ OptParams: This code is for analyzing and visualizing the L-curve and for 
        Frequency Analysis methods to find optimal reconstruction parameters 
        (alpha1 in TV based reconstructions). Please cite: 
        Milovic C, Prieto C, Bilgic B, Uribe S, Acosta-Cabronero J, Irarrazaval P, 
        Tejos C. Comparison of Parameter Optimization Methods for Quantitative 
        Susceptibility Mapping. Magn Reson Med. 2020. DOI: 10.1002/MRM.28435 

+ Relaxometry: This code is to estimate the initial magnetization and relaxation 
        time (T2/T2*) from MRI acquisitions. We provide a simple log fitting 
        funtion and nonlinear regularized functions. Regularization is done 
        explicitly on the relaxometry parameters, or implicitly on the projection 
        of the relaxometry parameters onto each simulated acquisition, for all
        echo times (see echoes suffix). These functions are designed for 3D data 
        acquired at multiple echo times only.

+ Unwrapping: We include both the analytic approximation and iterative methods 
        proposed by Schofield, M . Fast phase unwrapping algorithm for 
        interferometric applications. Opt Lett. 2003;28(14):1194-6 to perform 
        Laplacian-based unwrapping. We also include a Laplacian-based alternative, 
        where the gradient of the image is calculated with different phase offsets, 
        to create different wrapping artifacts, and a Poisson equation is solved 
        to find an optimal solution. Finally, an iterative implementation imposes 
        nonlinear data-fidelity on the acquired phase to avoid removing harmonic 
        components (similar to the PCG-Laplacian method, Robinson et al. 
        NMR Biomed 2016  DOI: 10.1002/nbm.3601). 
        Please refer to each file header for more details.
 


We provide an analytic Susceptibility brain phantom (data phantom folder) based 
on the phantom developed by C Langkammer, et al NeuroImage 2015. 
doi: 10.1016/j.neuroimage.2015.02.041, and C. Wisnieff, et al NeuroImage 2014 
doi: 10.1016/j.neuroimage.2012.12.050. 
Please cite our work (doi: 10.1002/mrm.27073) and theirs if you use this phantom.
 

Sample data from the 2016 QSM Reconstruction Challenge is also provided in a 
separate folder (data_challenge). In addition to the original single-orientation 
acquisition, COSMOS and X33 ground truths, we provide the X12 and X13 anisotropic 
components and a mask that estimates external non-tissue areas. If you use this 
dataset, please cite:
Langkammer, C., Schweser, F., Shmueli, K., Kames, C., Li, X., Guo, L., Milovic, C., 
Kim, J., Wei, H., Bredies, K., Buch, S., Guo, Y., Liu, Z., Meineke, J., Rauscher, A., 
Marques, J. P. and Bilgic, B. (2017), Quantitative susceptibility mapping: Report 
from the 2016 reconstruction challenge. Magn. Reson. Med. doi:10.1002/mrm.26830, 
and Milovic C, Tejos C, Acosta-Cabronero J, Ozbay PS, Schweser F, Marques JP, 
Irarrazaval P, Bilgic B, Langkammer C. The 2016 QSM Challenge: Lessons learned 
and considerations for a future challenge design; Magn Reson Med. 
2020;85:1624-1637 doi: 10.1002/mrm.28185

 
 
DISCLAIMER:
THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL CARLOS MILOVIC OR HIS 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY 
OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, BUSINESS INTERRUPTION; 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; AND LOSS OF USE, DATA OR PROFITS) 
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.