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
Usage Example:
      import scipy.io as load
      Signal = load.loadmat('SCGwithResp.mat')
      resp = Signal['resp'].ravel()
      from scaloGramDemo import scaloGramDemo
      Ta, Fa, St, PSD = scaloGramDemo(resp, 100, 0.1, 0.6, 0.005, True)
