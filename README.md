# PODforTumors
###  Functions and examples for reduced order modeling of PDEs applied to triple negative breast cancer
Associated manuscript: Fast model calibration for predicting the response of breast cancer to chemotherapy using proper orthogonal decomposition

Authors: Chase Christenson, Chengyue Wu, David A. Hormuth II, Graham Pash, Casey Stowers, Megan LaMonica, Karen Wilcox, Jingfei Ma, Gaiane M. Rauch, Thomas E. Yankeelov

Code authors: Chase Christenson, Graham Pash, Casey Stowers

The provided examples walk through the build procedure for a reduced order model and its subsequent usage to solve the forward model. This reduced forward model is used in a least squares minimization for parameter estimation to fit the model to the provided data examples.

The examples are set up to solve the following PDE (with variations):
$$\frac{\partial N}{\partial t}=\nabla \cdot \left(D(\textbf{x}) \nabla (N)\right)+k_p(\textbf{x})N\left(1-N\right) - \displaystyle\sum_{i=1}^{2}\alpha_i C(\textbf{x},t)e^{-\beta_i t} \tag{1}$$
For details on the model development:
C. Wu, A.M. Jarrett, Z. Zhou, et al., MRI-Based Digital Models Forecast Patient-Specific Treatment Responses to Neoadjuvant Chemotherapy in Triple-Negative Breast Cancer, Cancer Research 82(18) (2022) 3394–3404, https://doi.org/10.1158/0008-5472.CAN-22-1329.
