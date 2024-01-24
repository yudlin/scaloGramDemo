def scaloGramDemo(X, fs, flo, fup, fres, sFig=True):
    '''
    This program is used to derive the scalogram and time-averaged spectrum by
    complex Morlet wavelet.
    Designed by Yue-Der Lin on July 27, 2023.

    Input parameters:
        X  : The input signal to be analyzed (N by 1 array).
       fs  : The sampling frequency of signal (in Hz).
       flo  : Lower bound for spectrum analysis.
      fup  : Upper bound for spectrum analysis.
     fres  : Frequency resolution in the spectrum.
     sFig  : True = Show the figures, False = Do not show the figures
             (default = True).
    Output parameters:Outputs:
       Ta  : Times for spectrogram.
       Fa  : Frequencies for spectrogram and spectrum.
       St  : Spectrogram of X.
      PSD  : Averaged spectrum of X.
    Example:
      import scipy.io as load
      Signal = load.loadmat('SCGwithResp.mat')
      resp = Signal['resp'].ravel()
      from scaloGramDemo import scaloGramDemo
      Ta, Fa, St, PSD = scaloGramDemo(resp, 100, 0.1, 0.6, 0.005, True)
    '''
    # Import modules:
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np

    # Initialization:
    matplotlib.rcParams['font.family'] = 'Times'
    #   For complex Morlet wavelet:
    delta = 1/fs
    Fa = np.arange(flo, fup+fres, fres)  # Frequencies for analysis:
    S = np.divide(np.array([fs]), Fa)
    L = np.size(X)
    Ta = np.arange(0, L*delta, delta)

    # Complex Morlet wavelet transform:
    Coef = cMorl(X, S)

    St = np.abs((np.multiply(Coef, np.matrix.conjugate(Coef)))/L)
    #     PSD: Array of (lens(s), ).
    PSD = np.mean(St, axis=1)

    if sFig == True:
        # Spectrogram for X:
        plt.figure(figsize = (7.2, 5.4))
        plt.pcolormesh(Ta, Fa, St, cmap='ocean_r', shading='nearest')

        # To get current axis (gca) object.
        ax = plt.gca()
        #    width = 2 points.
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['top'].set_linewidth(2)
        ax.spines['left'].set_linewidth(2)
        ax.spines['right'].set_linewidth(2)

        plt.xticks(fontsize = 24, fontweight = 'bold')
        plt.yticks(fontsize = 24, fontweight = 'bold')
        plt.xlim(0, L/fs)
        plt.ylim(flo, fup)  # Show spectrogram in frequency range [flo, fup].
        plt.colorbar()
        plt.xlabel('Time (sec)', fontsize = 28, fontweight = 'bold')
        plt.ylabel('Frequency (Hz)', fontsize = 28, fontweight = 'bold')
        # plt.title('Spectrogram', fontsize = 15)
        plt.savefig('Fig_Spectrogram.jpg', dpi=1200, transparent=True,\
            bbox_inches='tight', pad_inches=0.1)

        # Averaged Spectrum for X:
        plt.figure(figsize = (7.2, 5.4))
        plt.plot(Fa, PSD, linewidth = 5)

        # To get current axis (gca) object.
        ax = plt.gca()
        #    width = 2 points.
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['top'].set_linewidth(2)
        ax.spines['left'].set_linewidth(2)
        ax.spines['right'].set_linewidth(2)

        plt.xticks(fontsize = 24, fontweight = 'bold')
        plt.yticks(fontsize = 24, fontweight = 'bold')
        plt.xlim(flo, fup)  # Show spectrum in frequency range [flo, fup].
        plt.xlabel('Frequency (Hz)', fontsize = 28, fontweight = 'bold')
        plt.ylabel('Power', fontsize = 28, fontweight = 'bold')
        # plt.title('Averaged Spectrum', fontsize = 15)
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))
        plt.savefig('Fig_Averaged_Spectrum.jpg', dpi=1200, transparent=True,\
            bbox_inches='tight', pad_inches=0.1)

        # Combination of Averaged Spectrum and Spectrogram for X:
        fig3 = plt.figure(figsize=(13, 6))
        plt.subplot(1,2,2)
        plt.pcolormesh(Ta, Fa, St, cmap = 'ocean_r', shading = 'nearest')
        plt.colorbar()

        # To get current axis (gca) object.
        ax = plt.gca()
        #    width = 2 points.
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['top'].set_linewidth(2)
        ax.spines['left'].set_linewidth(2)
        ax.spines['right'].set_linewidth(2)

        plt.xticks(fontsize = 20, fontweight = 'bold')
        plt.yticks(fontsize = 20, fontweight = 'bold')
        plt.xlabel('Time (sec)', fontsize = 24, fontweight = 'bold')
        # plt.ylabel(‘Frequency (Hz)', fontsize = 24, fontweight = 'bold')
        # plt.title('Spectrogram', fontsize = 15)
        plt.xlim(0, L/fs); plt.ylim(flo, fup)

        # fig3.subplots_adjust(wspace = 0.4)
        plt.subplot(1,2,1)
        plt.plot(PSD, Fa, linewidth = 3)
        # To get current axis (gca) object.
        ax = plt.gca()
        #    width = 2 points.
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['top'].set_linewidth(2)
        ax.spines['left'].set_linewidth(2)
        ax.spines['right'].set_linewidth(2)

        plt.xticks(fontsize = 20, fontweight = 'bold')
        plt.yticks(fontsize = 20, fontweight = 'bold')
        plt.xlabel('Power', fontsize = 24, fontweight = 'bold')
        plt.ylabel('Frequency (Hz)', fontsize = 24, fontweight = 'bold')
        # plt.title('Averaged Spectrum', fontsize = 15)
        plt.ylim(flo, fup)
        plt.ticklabel_format(axis='x', style='sci', scilimits=(-3,3))
        # fig3.subplots_adjust(wspace = 0.4)
        plt.savefig('Fig_Wavelet_Spectrum_Spectrogram.jpg', dpi=1200, transparent=True,\
            bbox_inches='tight', pad_inches=0.1)
    else:
        pass

    return Ta, Fa, St, PSD

def cMorl(Y, S):
    '''
    Complex Morlet wavelet transform: by Yue-Der Lin, on 2023/07/20.
    This function takes an input signal and computes the continuous wavelet
    transform at different scales using a complex Morlet Wavelet.

    References:
        C. Torrence and G. P. Compo, "A practical guide to wavelet analysis",
        Bull. Amer. Meteor. Soc., 1998, 79(1), pp.61-78.

    Inputs:
        Y  : The input signal to be analyzed (an array).
        S  : The scales coefficients for wavelet analysis (an array).
    Output:
     Coef  : The matrix returning CWT of input signal Y at each scale.
    '''
    # Import modules:
    import numpy as np
    import scipy
    from scipy import fft

    # Find power of 2 nearest to len(Y):
    N1 = len(Y)
    base2 = int(np.floor(np.log(N1)/np.log(2) + 0.4999))

    X = np.concatenate((Y, np.zeros(2**(base2+1)-N1)), axis=None)
    N = len(X)

    # Compute FFT of the padded time series:
    Fx = scipy.fft.fft(X)

    J1 = len(S) # Number of scales.
    k0 = 2*np.pi*1 # Center frequency of complex Morlet wavelet is 1.

    K = np.arange(np.floor(N/2))+1
    omega = 2*np.pi/N
    K = omega*K
    K1 = -K[int(np.floor((N-1)/2)):0:-1]
    K = np.concatenate((0., K), axis=None)
    # Eq.(5).
    K = np.concatenate((K, K1), axis=None)

    CoefT = np.zeros((J1, K.size), dtype=np.complex)
    for m in range(J1):
        expnt = -(S[m]*K-k0)**2/2*(K > 0)
        norm = np.sqrt(S[m]*omega)*(np.pi*(-0.25))*np.sqrt(N)
        Daughter = norm*np.exp(expnt)
        W = Daughter*(K > 0)
        CoefT[m,:] = scipy.fft.ifft(Fx*W) # Eq.(4).

    Coef =  CoefT[:, 0:N1]
    return Coef