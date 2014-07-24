
/**************************************************************************

                 Custom Added Filters

**************************************************************************/

const char *ITK_BINARY_DILATE_IMAGE_FILTER_DESC = "BinaryDilateImageFilter\n"
"Performs morphological dilation of the input binary image using a ball structuring element with a user-specified radius\n";

const char *ITK_BINARY_THINNING_IMAGE_FILTER_DESC="Binary Thinning Image Filter\n"
"The filter will produce a skeleton of the object.\n"
"The output background values are 0, and the foreground values are 1.\n"
"The input is assumed to be a binary image.\n";

const char *ITK_LAPLACIAN_RECURSIVE_GAUSSIAN_IMAGE_FILTER_DESC = "Laplacian of Gaussian (LoG) Filter\n"
"Computes the Laplacian of Gaussian (LoG) of an image by convolution with the second derivative of a Gaussian. This filter is implemented using the recursive gaussian filters.\n";

const char *ITK_SMOOTHING_RECURSIVE_GAUSSIAN_IMAGE_FILTER_DESC = "Smoothing Recursive Gaussian Image Filter\n"
"Computes the smoothing of an image by convolution with the Gaussian kernels implemented as IIR filters.\n"
"This filter is implemented using the recursive gaussian filters. For multi-component images, the filter works on each component independently.\n";

const char *ITK_INTERMODES_THRESHOLD_IMAGE_FILTER_DESC = "Prewitt's Intermodes Thresholdind Algorithm\n"
"\nAssumes a bimodal histogram. The histogram needs is smoothed (using a running average of size 3, iteratively) until there are only two local maxima. j and k Threshold t is (j+k)/2. Images with histograms having extremely unequal peaks or a broad and flat valley are unsuitable for this method.\n"
"\nJ. M. S. Prewitt and M. L. Mendelsohn, \"The analysis of cell images,\" in Annals of the New York Academy of Sciences, vol. 128, pp. 1035-1053, 1966\n";

const char *ITK_THRESHOLD_MAXIMUM_CONNECTED_COMPONENT_IMAGE_FILTER_DESC = "Pikaz's Topological Stable State Thresholding Algorithm\n"
"Implements a variant of Pikaz's Topological Stable State Thresholding Algorithm:\n"
"\n 1) Urish KL, August J, Huard J. \"Unsupervised segmentation for myofiber counting in immunoflourescent images\". Insight Journal. ISC/NA-MIC/MICCAI Workshop on Open-Source Software (2005)\n"
"\n 2)Pikaz A, Averbuch, A. \"Digital image thresholding based on topological stable-state\". Pattern Recognition, 29(5): 829-843, 1996.\n"
"The algorithm finds the threshold value that results in the maximum number of foregound connected components in the image that are larger than the specified minimal object size\n"
"SetMinimumObjectSizeInPixels() sets the minimum size for a connected component to be counted as an object\n";

const char *ITK_TRIANGLE_THRESHOLD_IMAGE_FILTER_DESC = "Traingle Threshold Algorithm\n"
"A line is drawn between the peak point in the hist and the furthest zero point (robustly estimated as the 1% or 99% point). The threshold is the position of maximum difference between the line and the original histogram.\n";

const char *ITK_RENYI_ENTROPY_THRESHOLD_IMAGE_FILTER_DESC = "Renyi Entropy-based Thresholding Algorithm\n"
"\nImplements Renyi's entropy based thresholding algorithm:\n"
"\nKapur J.N., Sahoo P.K., and Wong A.K.C. (1985) \"A New Method for Gray-Level Picture Thresholding Using the Entropy of the Histogram\" Graphical Models and Image Processing, 29(3): 273-285\n";

const char *ITKHUANGTHRESHOLDIMAGEFILTERDESC = "Huang's Fuzzy Similarity Thresholding Algorithm\n"
"\nImplements Huang's fuzzy similarity threshold algorithm:\n"
"\nHuang L.-K. and Wang M.-J.J. (1995) \"Image Thresholding by Minimizing the Measures of Fuzziness\" Pattern Recognition, 28(1): 41-51\n";

const char *ITKMOMENTSTHRESHOLDIMAGEFILTERDESC = "Tsai's Moment-preserving Thresholding Algorithm\n"
"\nImplements Tsai's moment preserving thresholding algorithm:\n"
"\nW. Tsai, \"Moment-preserving thresholding: a new approach,\" Computer Vision, Graphics, and Image Processing, vol. 29, pp. 377-393, 1985.\n";

const char *ITKISODATATHRESHOLDIMAGEFILTERDESC = "IsoData Threshold Algorithm\n"
"\nImplements Ridler and Calvard's IsoData Threshold Algorithm:\n"
"T.W. Ridler, S. Calvard, Picture thresholding using an iterative selection method, IEEE Trans. System, Man and Cybernetics, SMC-8 (1978) 630-632\n"
"The procedure divides the image into objects and background by taking an initial threshold, then the averages of the pixels at or below the threshold and pixels above are computed. The averages of those two values are computed, the threshold is incremented and the process is repeated until the threshold is larger than the composite average. That is, threshold = (average background + average objects)/2\n";

const char *ITKYENTHRESHOLDIMAGEFILTERDESC = "Yen's Entropy Based Thresholding Algorithm\n"
"\nImplements Yen thresholding method\n"
"1) Yen J.C., Chang F.J., and Chang S. (1995) \"A New Criterion for Automatic Multilevel Thresholding\" IEEE Trans. on Image Processing, 4(3): 370-378\n"
"2) Sezgin M. and Sankur B. (2004) \"Survey over Image Thresholding Techniques and Quantitative Performance Evaluation\" Journal of Electronic Imaging, 13(1): 146-165\n";

const char *ITKSHANBHAGTHRESHOLDIMAGEFILTERDESC = "Shanbag's Fuzzy Entropic Thresholding\n"
"\nImplements the Shanbhag's Fuzzy Entropic Thresholding algorithm\n"
"Shanhbag A.G. (1994) \"Utilization of Information Measure as a Means of Image Thresholding\" Graphical Models and Image Processing, 56(5): 414-419\n";

const char *ITKLITHRESHOLDIMAGEFILTERDESC = "Li's Maximum Cross Entropy Threshold\n"
"\nImplements Li's Minimum Cross Entropy thresholding method This implementation is based on the iterative version (Ref. 2) of the algorithm.\n"
"1) Li C.H. and Lee C.K. (1993) \"Minimum Cross Entropy Thresholding\" Pattern Recognition, 26(4): 617-625\n" 
"2) Li C.H. and Tam P.K.S. (1998) \"An Iterative Algorithm for Minimum Cross Entropy Thresholding\" Pattern Recognition Letters, 18(8): 771-776\n" 
"3) Sezgin M. and Sankur B. (2004) \"Survey over Image Thresholding Techniques and Quantitative Performance Evaluation\" Journal of Electronic Imaging, 13(1): 146-165\n";

const char *ITKMAXIMUMENTROPYTHRESHOLDIMAGEFILTERDESC = "Maximum Entropy Threshold\n"
"Thresholds the image using Kapur and Sahoo's Maximum Entropy Threshold Algorithm\n"
"Implements Kapur-Sahoo-Wong (Maximum Entropy) thresholding method Kapur J.N., Sahoo P.K., and Wong A.K.C. (1985) \"A New Method for Gray-Level Picture Thresholding Using the Entropy of the Histogram\" Graphical Models and Image Processing, 29(3): 273-285\n";

const char *ITKOTSUTHRESHOLDIMAGEFILTERDESC2 = "Otsu's Thresholding Algorithm\n"
"Thresholds an image using the Otsu's Thresholding algorithm\n"
"This filter creates a binary thresholded image that separates an image into foreground and background components. The filter computes the threshold using the OtsuThresholdCalculator and applies that theshold to the input image using the BinaryThresholdImageFilter.";

const char *ITKKITTLERILLINGWORTHTHRESHOLDIMAGEFILTERDESC = "Minimum Error Thresholding\n"
"Thresholds the image using Kittler and Illingworth's Minimum Error Thresholding Algorithm.\n"
"Kittler and J. Illingworth, \"Minimum error thresholding,\" Pattern Recognition, vol. 19, pp. 41-47, 1986. C. A. Glasbey, \"An analysis of histogram-based thresholding algorithms,\" CVGIP: Graphical Models and Image Processing, vol. 55, pp. 532-537, 1993";

const char *ITKGRADIENTANISOTROPICDIFFUSIONIMAGEFILTERDESC = "Gradient Anisotropic Diffusion Image Filter\n"
"This filter performs anisotropic diffusion on a scalar itk::Image using the classic Perona-Malik, gradient magnitude based equation implemented in itkGradientNDAnisotropicDiffusionFunction\n";

const char *ITKMEDIANIMAGEFILTERDESC = "Median Image Filter\n"
" Applies a median filter to an image.\n"
"\n"
" Computes an image where a given pixel is the median value of the the pixels in a neighborhood about the corresponding input pixel.\n"
"\n"
" A median filter is one of the family of nonlinear filters. It is used to smooth an image without being biased by outliers or shot noise.\n"
"\n"
"SetRadius specifies the neighborhood size of the median filter";

const char *ITKOTSUTHRESHOLDIMAGEFILTERDESC = "Otsu Threshold Image Filter";

const char *ITKANTIALIASBINARYIMAGEFILTERDESC = "Anti-Alias Binary Image Filter";

const char *ITKSIGNEDMAURERDISTANCEMAPIMAGEFILTERDESC = "Signed Maurer Distance Map Image Filter";

const char *ITKPOINTSETVOLREGISTRATIONDESC = "Registration of two volumes based on point sets";

const char *ITKLABELEDVOLREGISTRATIONDESC = "Registration of two labeled volumes"; 

//////////////////Description from itkDerivativeImageFilter.h//////////////////
const char *ITKDERIVATIVEIMAGEFILTERDESC=" DerivativeImageFilter\n"
"  Computes the directional derivative of an image.\n"
" The directional derivative at each pixel location is computed by convolution\n"
" with a derivative operator of user-specified order.\n"
"\n"
" SetOrder specifies the order of the derivative.\n"
"\n"
" SetDirection specifies the direction of the derivative with respect to the\n"
" coordinate axes of the image.\n"
" \n"
"  Image\n"
"  Neighborhood\n"
"  NeighborhoodOperator\n"
"  NeighborhoodIterator\n"
" \n"
"  ImageFeatureExtraction";


//////////////////Description from itkDemonsRegistrationFilter.h//////////////////
const char *ITKDEMONSREGISTRATIONFILTERDESC=" DemonsRegistrationFilter\n"
"  Deformably register two images using the demons algorithm.\n"
"\n"
" DemonsRegistrationFilter implements the demons deformable algorithm that \n"
" register two images by computing the deformation field which will map a \n"
" moving image onto a fixed image.\n"
"\n"
" A deformation field is represented as a image whose pixel type is some\n"
" vector type with at least N elements, where N is the dimension of\n"
" the fixed image. The vector type must support element access via operator\n"
" []. It is assumed that the vector elements behave like floating point\n"
" scalars.\n"
"\n"
" This class is templated over the fixed image type, moving image type\n"
" and the deformation field type.\n"
"\n"
" The input fixed and moving images are set via methods SetFixedImage\n"
" and SetMovingImage respectively. An initial deformation field maybe set via\n"
" SetInitialDeformationField or SetInput. If no initial field is set,\n"
" a zero field is used as the initial condition.\n"
"\n"
" The algorithm has one parameters: the number of iteration to be performed.\n"
"\n"
" The output deformation field can be obtained via methods GetOutput\n"
" or GetDeformationField.\n"
"\n"
" This class make use of the finite difference solver hierarchy. Update\n"
" for each iteration is computed in DemonsRegistrationFunction.\n"
"\n"
"  This filter assumes that the fixed image type, moving image type\n"
" and deformation field type all have the same number of dimensions.\n"
" \n"
"  DemonsRegistrationFunction \n"
"  DeformableImageRegistration MultiThreaded";

//////////////////Description from itkThinPlateSplineKernelTransform.h//////////////////
const char *ITKTHINPLATESPLINEKERNELTRANSFORMDESC=" ThinPlateSplineKernelTransform\n"
" This class defines the thin plate spline (TPS) transformation.\n"
" It is implemented in as straightforward a manner as possible from\n"
" the IEEE TMI paper by Davis, Khotanzad, Flamig, and Harms,\n"
" Vol. 16 No. 3 June 1997\n"
"\n"
"  Transforms";

//////////////////Description from itkGradientMagnitudeImageFilter.h//////////////////
const char *ITKGRADIENTMAGNITUDEIMAGEFILTERDESC=" GradientMagnitudeImageFilter\n"
"  Computes the gradient magnitude of an image region at each pixel.\n";

//////////////////Description from itkGradientMagnitudeRecursiveGaussianImageFilter.h//////////////////
const char *ITKGRADIENTMAGNITUDERECURSIVEGAUSSIANIMAGEFILTERDESC=" GradientMagnitudeRecursiveGaussianImageFilter\n"
"  Computes the Magnitude of the Gradient of an image by convolution\n"
"        with the first derivative of a Gaussian.\n";

