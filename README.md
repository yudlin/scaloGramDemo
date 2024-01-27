For scaloGramDemo.py:

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
      scg = Signal['scg'].ravel()
      from scipy import signal
      resp = signal.resample(resp, int(len(resp)*5/100))
      scg = signal.resample(scg, int(len(scg)*5/100))
      from scaloGramDemo import scaloGramDemo
      Ta1, Fa1, St1, PSD1 = scaloGramDemo(resp, 5, 0.05, 0.6, 0.005, True)
      Ta2, Fa2, St2, PSD2 = scaloGramDemo(scg, 5, 0.05, 0.6, 0.005, True)


For Test Datsets:

SCGwithResp.mat: Test data for SCG and respiration signal (from CEBS database).

f04sSCGwithResp.mat: Test data for SCG and respiration signal (from experiments in FCU).

fingerPPGwithRIIV.mat: Test data for finger PPG and respiration signal (from MIMIC II database).

wristPPGwithRIIV.mat: Test data for wrist PPG and respiration signal (from experiments in FCU). 
