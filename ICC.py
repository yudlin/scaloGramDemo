def ICC(M, type, alpha=0.05, r0=0):
    '''
    This code is converted from Arash Salarian's MATLAB code, which is used to estimate the
    intraclass correlation (ICC) between two arrays.

    M is matrix of observations. Each row is an object of measurement and each column is a
    judge or measurement.

    'type' is a string that can be one of the six possible codes for the
    desired type of ICC:
      '1-1': The degree of absolute agreement among measurements made on randomly selected
             objects. It estimates the correlation of any two measurements.
      '1-k': The degree of absolute agreement of measurements that are averages of k
             independent measurements on randomly selected objects.
      'C-1': case 2: The degree of consistency among measurements. Also known as norm-referenced
             reliability and as Winer's adjustment for anchor points.
             case 3: The degree of consistency among measurements made under the fixed levels of
             the column factor.
             This ICC estimates the corrlation of any two measurements, but when interaction is
             present, it underestimates reliability.
      'C-k': case 2: The degree of consistency for measurements that are averages of k independent
             measurements on randomly selected on bjectgs. Known as Cronbach's alpha in psychometrics.
             case 3: The degree of consistency for averages of k independent measures made under the
             fixed levels of column factor.
      'A-1': case 2: The degree of absolute agreement among measurements. Also known as criterion-
             referenced reliability.
             case 3: The absolute agreement of measurements made under the fixed levels of the column
             factor.
      'A-k': case 2: The degree of absolute agreement for measurements that are averages of k
             independent measurements on randomly selected objects.
             case 3: he degree of absolute agreement for measurements that are based on k independent
             measurements maded under the fixed levels of the column factor.

    Reference: McGraw, K. O., Wong, S. P., "Forming Inferences About Some Intraclass
        Correlation Coefficients", Psychological Methods, Vol. 1, No. 1, pp. 30-46, 1996.

    Example usage:
### Replace 'M' with your actual data matrix, and choose the appropriate 'type'
import numpy as np
from numpy import random
from ICC import ICC
R1 = np.random.rand(10,1); R2 = np.random.rand(10,1)
M = np.column_stack((R1, R2));
type = 'A-1'
r, LB, UB, F, df1, df2, p = ICC(M, type)
# r, LB, UB, _, _, _, _ = ICC(M, type)
print(f"Intraclass correlation (ICC): {r}")
print(f"Lower bound (LB): {LB}")
print(f"Upper bound (UB): {UB}")
print(f"F-value: {F}")
print(f"Degree of freedom 1: {df1}")
print(f"Degree of freedom 2: {df2}")
print(f"P-value: {p}")
    '''
    import numpy as np
    from scipy.stats import f

    n, k = M.shape

    # SStotal = np.var(M.flatten()) * (n*k - 1)
    SStotal = np.var(M.flatten())*(n*k - 1)
    MSR = np.var(np.mean(M, axis=1)) * k
    MSW = np.sum(np.var(M, axis=1, ddof=0))/n
    MSC = np.var(np.mean(M, axis=0))*n
    MSE = (SStotal - MSR*(n - 1) - MSC*(k - 1))/((n - 1)*(k - 1))

    def ICC_case_1_1(MSR, MSE, MSC, MSW, alpha, r0, n, k):
        r = (MSR - MSW)/(MSR + (k-1)*MSW)

        F = (MSR/MSW)*(1-r0)/(1+(k-1)*r0)
        df1 = n-1
        df2 = n*(k-1)
        p = 1 - f.cdf(F, df1, df2)

        FL = (MSR/MSW)/f.ppf(1-alpha/2, df1, df2)
        FU = (MSR/MSW)*f.ppf(1-alpha/2, df2, df1)

        LB = (FL - 1)/(FL + (k-1))
        UB = (FU - 1)/(FU + (k-1))

        return r, LB, UB, F, df1, df2, p

    def ICC_case_1_k(MSR, MSE, MSC, MSW, alpha, r0, n, k):
        r = (MSR - MSW)/MSR

        F = (MSR/MSW)*(1-r0)
        df1 = n-1
        df2 = n*(k-1)
        p = 1 - f.cdf(F, df1, df2)

        FL = (MSR/MSW)/f.ppf(1-alpha/2, df1, df2)
        FU = (MSR/MSW)*f.ppf(1-alpha/2, df2, df1)

        LB = 1 - 1/FL
        UB = 1 - 1/FU

        return r, LB, UB, F, df1, df2, p

    def ICC_case_C_1(MSR, MSE, MSC, MSW, alpha, r0, n, k):
        r = (MSR - MSE)/(MSR + (k-1)*MSE)

        F = (MSR/MSE)*(1-r0)/(1+(k-1)*r0)
        df1 = n - 1
        df2 = (n-1)*(k-1)
        p = 1 - f.cdf(F, df1, df2)

        FL = (MSR/MSE)/f.ppf(1-alpha/2, df1, df2)
        FU = (MSR/MSE)*f.ppf(1-alpha/2, df2, df1)

        LB = (FL - 1)/(FL + (k-1))
        UB = (FU - 1)/(FU + (k-1))

        return r, LB, UB, F, df1, df2, p

    def ICC_case_C_k(MSR, MSE, MSC, MSW, alpha, r0, n, k):
        r = (MSR - MSE)/MSR

        F = (MSR/MSE)*(1-r0)
        df1 = n - 1
        df2 = (n-1)*(k-1)
        p = 1 - f.cdf(F, df1, df2)

        FL = (MSR/MSE)/f.ppf(1-alpha/2, df1, df2)
        FU = (MSR/MSE)*f.ppf(1-alpha/2, df2, df1)

        LB = 1 - 1/FL
        UB = 1 - 1/FU

        return r, LB, UB, F, df1, df2, p

    def ICC_case_A_1(MSR, MSE, MSC, MSW, alpha, r0, n, k):
        r = (MSR - MSE)/(MSR + (k-1)*MSE + k*(MSC-MSE)/n)

        a = (k*r0)/(n*(1-r0))
        b = 1 + (k*r0*(n-1))/(n*(1-r0))
        F = MSR/(a*MSC + b*MSE)

        a = k*r/(n*(1-r))
        b = 1+k*r*(n-1)/(n*(1-r))
        v = (a*MSC + b*MSE)**2/((a*MSC)**2/(k-1) + (b*MSE)**2/((n-1)*(k-1)))

        df1 = n - 1
        df2 = v
        p = 1 - f.cdf(F, df1, df2)

        Fs = f.ppf(1-alpha/2, df1, df2)
        LB = n*(MSR - Fs*MSE)/(Fs*(k*MSC + (k*n - k - n)*MSE) + n*MSR)

        Fs = f.ppf(1-alpha/2, df2, df1)
        UB = n*(Fs*MSR-MSE)/(k*MSC + (k*n - k - n)*MSE + n*Fs*MSR)

        return r, LB, UB, F, df1, df2, p

    def ICC_case_A_k(MSR, MSE, MSC, MSW, alpha, r0, n, k):
        r = (MSR - MSE)/(MSR + (MSC-MSE)/n)

        c = r0/(n*(1-r0))
        d = 1 + (r0*(n-1))/(n*(1-r0))
        F = MSR/(c*MSC + d*MSE)

        a = k*r/(n*(1-r))
        b = 1+k*r*(n-1)/(n*(1-r))
        v = (a*MSC + b*MSE)**2/((a*MSC)**2/(k-1) + (b*MSE)**2/((n-1)*(k-1)))

        df1 = n - 1
        df2 = v
        p = 1 - f.cdf(F, df1, df2)

        Fs = f.ppf(1-alpha/2, df1, df2)
        LB = n*(MSR - Fs*MSE)/(Fs*(MSC-MSE) + n*MSR)

        Fs = f.ppf(1-alpha/2, df2, df1)
        UB = n*(Fs*MSR - MSE)/(MSC - MSE + n*Fs*MSR)

        return r, LB, UB, F, df1, df2, p

    switch = {
        '1-1': ICC_case_1_1,
        '1-k': ICC_case_1_k,
        'C-1': ICC_case_C_1,
        'C-k': ICC_case_C_k,
        'A-1': ICC_case_A_1,
        'A-k': ICC_case_A_k,
    }

    if type in switch:
        return switch[type](MSR, MSE, MSC, MSW, alpha, r0, n, k)
    else:
        raise ValueError(f"Invalid type: {type}. Supported types are '1-1', '1-k', 'C-1', 'C-k', 'A-1', 'A-k'.")
