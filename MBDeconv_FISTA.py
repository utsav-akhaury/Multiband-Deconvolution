import numpy as np
from scipy.signal import convolve
from skimage.transform import rescale
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

def nmse(signal_1, singal_2):

    r"""
    Compute the Normalised Mean Square Error (NMSE) between two images.

    Parameters
    ----------
    signal_1 : 2D image 1
    signal_2 : 2D image 2

    Returns
    ------- 
    nmse_val : float
        The NMSE between the two images.
    """
    nmse_val = (np.linalg.norm(singal_2 - signal_1, axis=(0,1))**2 /
                np.linalg.norm(signal_1, axis=(0,1))**2)
    return nmse_val


def max_sv(lsst_psf, euc_psf, sigma_noise, lamdbd_constr, SED):

    r"""
    Compute the square of the spectral radius of the convolution matrix.

    Parameters
    ----------
    lsst_psf : 2D array
        The LSST PSF at Euclid resolution.
    euc_psf : 2D array
        The Euclid PSF.
    sigma_noise : 1D array with 4 float values (3 for LSST bands and 1 for Euclid VIS)
        The noise level in each band.
    lamdbd_constr : float   
        The constraint parameter.
    SED : 1D array with 3 float values (for each LSST band)
        The fractional contribution of each LSST band to the Euclid VIS band.

    Returns
    -------
    sv : 1D array with 4 float values (3 for LSST bands and 1 for Euclid VIS)
        The square of the spectral radius of the convolution matrix for each band.
    """

    H_lsst_psf = np.fft.fft2(lsst_psf, axes=(0,1))
    normH_lsst = np.abs(np.rot90(H_lsst_psf, 2, axes=(0,1))*H_lsst_psf) 

    H_psf_euc = np.fft.fft2(euc_psf)
    normH_euc = np.abs(np.rot90(H_psf_euc, 2)*H_psf_euc)
    
    res = (normH_lsst / np.linalg.norm(sigma_noise[:-1])**2 +
           2. * lamdbd_constr * SED**2 * np.expand_dims(normH_euc, axis=-1) / np.linalg.norm(sigma_noise[-1])**2)
    
    sv = np.real(np.amax(res, axis=(0,1)))
    return sv


def get_alpha(sv):

    r"""
    Compute the gradient step size.

    Parameters
    ----------
    sv : float
        The square of the spectral radius of the convolution matrix.

    Returns
    -------
    alpha : 1D array
        The gradient step size.
    """

    alpha = 1.0 / (sv * (1.0 + 1.0e-5))
    return alpha


def proj(x):

    r"""
    Enforce non-negative values.

    Parameters
    ----------
    x : N-D array
        The input array.

    Returns
    -------
    The output array with non-negative values. : N-D array
    """

    return np.maximum(np.real(x), 0)


def rescale_preserve_flux(data, rescale_fact):

    r"""
    Rescale an image while preserving the flux.

    Parameters
    ----------
    data : 2D array
        The input image.

    rescale_fact : float
        The rescale factor.

    Returns
    -------
    The rescaled image : 2D array
    """

    if rescale_fact==1:
        return data 
    else:
        arr_flux = np.sum(data)
        arr_rescaled = rescale(data, rescale_fact, mode='constant', anti_aliasing=True)
        return (arr_rescaled / np.sum(arr_rescaled)) * arr_flux
    

def H(data, psf, rescale_fact):

    r"""
    Convolve an image with a PSF (H operation).
        
    Parameters
    ----------
    data : 2D array
        The input image.
    psf : 2D array
        The PSF.
    rescale_fact : float
        The rescale factor.
            
    Returns
    -------
    The convolved image : 2D array
    """

    return rescale_preserve_flux(convolve(data, psf, mode='same'), rescale_fact)


def Ht(data, psf, rescale_fact):

    r"""
    Convolve an image with a PSF rotated by 180_degrees (H_transpose operation).

    Parameters
    ----------
    data : 2D array
        The input image.
    psf : 2D array
        The PSF.
    rescale_fact : float
        The rescale factor.

    Returns
    -------
    The convolved image : 2D array
    """
    
    return rescale_preserve_flux(convolve(data, np.rot90(psf, 2), mode='same'), rescale_fact)


def grad(y_lsst, y_euc, ch, x_k, lsst_psf, euc_psf, sigma_noise, lamdbd_constr, SED):

    r"""
    Compute the gradient of the loss function.
    
    Parameters
    ----------
    y_lsst : 3D array
        The LSST image in each LSST band.
    y_euc : 2D array
        The Euclid image in the Euclid VIS band.
    ch : int
        The band index.
    x_k : 3D array
        The current estimate of the deconvolved image in each LSST band.
    lsst_psf : 3D array
        The LSST PSF in each LSST band at Euclid resolution.
    euc_psf : 2D array  
        The Euclid PSF.
    sigma_noise : 1D array with 4 float values (3 for LSST bands and 1 for Euclid VIS)
        The noise level in each band.
    lamdbd_constr : float
        The constraint parameter.
    SED : 1D array with 3 float values (for each LSST band)
        The fractional contribution of each LSST band to the Euclid VIS band.
        
    Returns
    -------
    The gradient of the loss function : 3D array
    """
    
    dat_fid = Ht(H(x_k[...,ch], lsst_psf[...,ch], 0.5) - y_lsst[...,ch], lsst_psf[...,ch], 2) / np.linalg.norm(sigma_noise[ch])**2

    constr_term = H(x_k[...,0]*SED[0] + 
                    x_k[...,1]*SED[1] + 
                    x_k[...,2]*SED[2], euc_psf, 1) - y_euc

    constr_grad = 2 * lamdbd_constr * SED[ch] * Ht(constr_term, euc_psf, 1) / np.linalg.norm(sigma_noise[-1])**2

    return dat_fid + constr_grad


def loss_func(y_lsst, y_euc, x_k, lsst_psf, euc_psf, sigma_noise, lamdbd_constr, SED):

    r"""
    Compute the loss function.

    Parameters
    ----------
    y_lsst : 3D array
        The LSST image in each LSST band.
    y_euc : 2D array
        The Euclid image in the Euclid VIS band.
    x_k : 3D array
        The current estimate of the deconvolved image in each LSST band.
    lsst_psf : 3D array
        The LSST PSF in each LSST band at Euclid resolution.
    euc_psf : 2D array
        The Euclid PSF.
    sigma_noise : 1D array with 4 float values (3 for LSST bands and 1 for Euclid VIS)
        The noise level in each band.
    lamdbd_constr : float
        The constraint parameter.
    SED : 1D array with 3 float values (for each LSST band)
        The fractional contribution of each LSST band to the Euclid VIS band.  

    Returns
    -------
    The loss function : 1D array
    """
    
    dat_fid = np.zeros((y_lsst.shape[-1]))
    dat_fid[0] = 0.5 * np.linalg.norm((H(x_k[...,0], lsst_psf[...,0], 0.5) - y_lsst[...,0]) / sigma_noise[0])**2
    dat_fid[1] = 0.5 * np.linalg.norm((H(x_k[...,1], lsst_psf[...,1], 0.5) - y_lsst[...,1]) / sigma_noise[1])**2
    dat_fid[2] = 0.5 * np.linalg.norm((H(x_k[...,2], lsst_psf[...,2], 0.5) - y_lsst[...,2]) / sigma_noise[2])**2
   
    constr = lamdbd_constr * np.linalg.norm((y_euc - H(np.dot(x_k, SED), euc_psf, 1)) / sigma_noise[-1], axis=(0,1))**2

    return dat_fid + constr


def set_lambda_constr(y_lsst, y_euc, x_k, lsst_psf, euc_psf, sigma_noise, SED, frac_contr=1.):

    r"""
    Compute the lambda constraint parameter (the multiplicative factor that gives weight to the constraint term).

    Parameters
    ----------
    y_lsst : 3D array
        The LSST image in each LSST band.
    y_euc : 2D array
        The Euclid image in the Euclid VIS band.
    x_k : 3D array
        The current estimate of the deconvolved image in each LSST band.
    lsst_psf : 3D array
        The LSST PSF in each LSST band at Euclid resolution.
    euc_psf : 2D array
        The Euclid PSF.
    sigma_noise : 1D array with 4 float values (3 for LSST bands and 1 for Euclid VIS)
        The noise level in each band.
    SED : 1D array with 3 float values (for each LSST band)
        The fractional contribution of each LSST band to the Euclid VIS band.
    frac_contr : float  
        The constraint parameter.

    Returns
    -------
    lamdbd_constr : 1D array
        The lambda constraint parameter.
    """


    dat_fid = np.zeros((y_lsst.shape[-1]))
    dat_fid[0] = 0.5 * np.linalg.norm((H(x_k[...,0], lsst_psf[...,0], 0.5) - y_lsst[...,0]) / sigma_noise[0])**2
    dat_fid[1] = 0.5 * np.linalg.norm((H(x_k[...,1], lsst_psf[...,1], 0.5) - y_lsst[...,1]) / sigma_noise[1])**2
    dat_fid[2] = 0.5 * np.linalg.norm((H(x_k[...,2], lsst_psf[...,2], 0.5) - y_lsst[...,2]) / sigma_noise[2])**2

    lamdbd_constr = frac_contr * (dat_fid / np.linalg.norm((y_euc - H(np.dot(x_k, SED), euc_psf, 1)) / sigma_noise[-1], axis=(0,1))**2) 

    return lamdbd_constr


def runDeconv(y_lsst, y_euc, lsst_psf, euc_psf, sigma_noise, n_iter, SED, frac_contr, lamdbd_constr=0):

    r"""
    Run the Multi-Channel Deconvolution algorithm.

    Parameters
    ----------
    y_lsst : 3D array
        The LSST image in each LSST band.
    y_euc : 2D array
        The Euclid image in the Euclid VIS band.
    lsst_psf : 3D array
        The LSST PSF in each LSST band at Euclid resolution.
    euc_psf : 2D array
        The Euclid PSF.  
    sigma_noise : 1D array with 4 float values (3 for LSST bands and 1 for Euclid VIS)
        The noise level in each band.
    n_iter : int
        The number of iterations.
    SED : 1D array with 4 float values (3 for LSST bands and 1 for Euclid VIS)
        The fractional contribution of each LSST band to the Euclid VIS band.
    frac_contr : float
        The constraint parameter.
    lamdbd_constr : float
        The constraint parameter.

    Returns
    -------
    x_k : 3D array
        The deconvolved image in each LSST band.
    """
    
    # Initialize the first estimate of the deconvolved images as the LSST images at Euclid resolution
    x_0 = np.zeros((y_euc.shape[0], y_euc.shape[1], y_lsst.shape[-1]))
    for j in range(y_lsst.shape[-1]):
        x_0[:,:,j] = rescale_preserve_flux(y_lsst[:,:,j], 2)

    # FISTA parameters
    x_k = x_0
    t_k = 1.0
          
    # square of spectral radius of convolution matrix
    sv = max_sv(lsst_psf, euc_psf, sigma_noise, lamdbd_constr, SED)           

    # The gradient descent step
    alpha = get_alpha(sv)/10

    for k in range(n_iter):

        # estimate lambd_constr
        lamdbd_constr = set_lambda_constr(y_lsst, y_euc, x_k, lsst_psf, euc_psf, sigma_noise, SED, frac_contr)
        
        ## Gradient Descent update
        x_k[...,0] = x_k[...,0] - alpha[0] * grad(y_lsst, y_euc, 0, x_k, lsst_psf, euc_psf, sigma_noise, lamdbd_constr[0], SED)
        x_k[...,1] = x_k[...,1] - alpha[1] * grad(y_lsst, y_euc, 1, x_k, lsst_psf, euc_psf, sigma_noise, lamdbd_constr[1], SED)
        x_k[...,2] = x_k[...,2] - alpha[2] * grad(y_lsst, y_euc, 2, x_k, lsst_psf, euc_psf, sigma_noise, lamdbd_constr[2], SED)
                
        # FISTA Update
        x_k1 = x_k
        x_k1 = proj(x_k1)
        t_k1 = (1. + np.sqrt(4.*t_k**2 + 1.))/2.
        lambda_fista = 1 + (t_k - 1)/t_k1
        x_k1 = x_k + lambda_fista*(x_k1 - x_k)  
        x_k = x_k1

    return x_k


def comparison_normcbar(targets, noisy, euclid, deconv, labels, figsize):

    r"""
    Plot the comparison between the noisy, deconvolved, and target images.
    
    Parameters
    ----------
    targets : 3D array
        The target images.
    noisy : 3D array
        The noisy images.
    euclid : 2D array
        The Euclid image.
    deconv : 3D array
        The deconvolved images.
    labels : list of strings
        The labels for the images.
    figsize : tuple
        The figure size.

    Returns
    -------
    fig : figure
        The figure with the comparison between the noisy, deconvolved, and target images.
    """
    
    list_im = []

    n_row = targets.shape[-1]
    n_col = len(labels)
    
    bands = ['$r$','$i$','$z$']
    bands_hst = ['$F606W$','$F775W$','$F850LP$']
    
    for i in range(n_row):       
        list_im += [noisy[...,i], 
                    deconv[...,i], 
                    targets[...,i]]

    fig = plt.figure(figsize=figsize)
    gs1 = gridspec.GridSpec(n_row, n_col+1)
    gs1.update(wspace=0.01, hspace=0.01)
        
    vmin = np.inf
    vmax = 0

    for t in range(n_row):
        vmin = np.min((vmin, np.min(list_im[t*n_col + 1])))
        vmax = np.max((vmax, np.max(list_im[t*n_col + 1])))

    for i in range(n_row):
        for k in range(n_col): 
            
            vmin = np.min((vmin, np.min(list_im[i*n_col + 1])))
            vmax = np.max((vmax, np.max(list_im[i*n_col + 1])))

            axes = plt.subplot(gs1[i*n_col + k + 1*(i+1)])
            axes.axis('off')

            if k==1:
                im = axes.imshow(list_im[i*n_col + k], origin='lower', vmin=vmin, vmax=vmax, cmap='afmhot')
                axes.text(.01, .938, labels[k]+': '+bands[i]+'-band', ha='left', va='bottom', fontsize=25, fontweight="bold", color='black', bbox={'facecolor': 'white', 'alpha': 1, 'pad': 5}, transform=axes.transAxes)
            elif k==2:
                im = axes.imshow(list_im[i*n_col + k], origin='lower', vmin=vmin, vmax=vmax, cmap='afmhot')
                axes.text(.01, .943, labels[k]+': '+bands_hst[i], ha='left', va='bottom', fontsize=25, fontweight="bold", color='black', bbox={'facecolor': 'white', 'alpha': 1, 'pad': 5}, transform=axes.transAxes)
            elif k==0:
                im = axes.imshow(list_im[i*n_col + k]/4, origin='lower', vmin=vmin, vmax=vmax, cmap='afmhot')
                axes.text(.01, .943, labels[k]+': '+bands[i]+'-band', ha='left', va='bottom', fontsize=25, fontweight="bold", color='black', bbox={'facecolor': 'white', 'alpha': 1, 'pad': 5}, transform=axes.transAxes)

    axes = plt.subplot(gs1[4])
    im = axes.imshow(euclid, origin='lower',  vmax=vmax, cmap='afmhot')
    axes.text(.01, .943, 'Euclid: $VIS$', ha='left', va='bottom', fontsize=25, fontweight="bold", color='black', bbox={'facecolor': 'white', 'alpha': 1, 'pad': 5}, transform=axes.transAxes) 
    axes.axis('off')

    plt.tight_layout()
    
    return fig