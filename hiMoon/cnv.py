"""
Copy number variation algorithm. Adapted from the DOC method using maximum penalized
likelihood estimation (MPLE) described by Chen, et al. 

https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1566-3#MOESM1

All functions are jitted with numba/LLVM. Should all run multithreaded
and fairly fast.
"""

import numpy as np

from numba import int32, float32, jit, void


@jit(float32(int32), nopython=True, parallel=True)
def calc_penalty(n: int) -> float:
    """
    Penalization factor (lambda)
    n is total reads in entire segment
    Based on Bayesian information criterion (BIC)
    TODO: Add AIC, P value, others?
    """
    return 0.5 * np.log(n)


@jit(float32(int32, int32, float32), nopython=True, parallel=True)
def log_lik(t: int,c: int,l: float) -> float:
    """
    Log likelihood calculation
    """
    if c*t != 0: # Make sure that neither c or t is 0, otherwise return 0
        p = (l * t) / (l * t + c) # Probability that a given read comes from the case
        return  l * t * np.log(p) + c * np.log(1 - p) # log-likelihood for the observed bernoulli trial
    else:
        return 0


@jit(int32(int32, int32, int32[:]), nopython=True, parallel=True)
def grc_on(j_pos: int, i_pos: int, reads: int) -> int:
    """
    j_pos = possible CNV start
    i_pos = possible CNV end
    reads = either case or control numpy array
    Returns the number of reads containing the CBP
    """
    return ((j_pos <= reads) & (i_pos >= reads)).sum()


@jit(int32(int32, int32, int32[:]), nopython=True, parallel=True)
def grc_off(j_pos: int, i_pos: int, reads: int) -> int:
    """
    j_pos = possible CNV start
    i_pos = possible CNV end
    reads = either case or control numpy array
    Returns the number of reads not containing the CBP
    """
    return ((j_pos > reads) | (i_pos < reads)).sum()


@jit(float32[:,:](int32[:], int32[:], int32), nopython=True, parallel=True)
def find_cnv(case: np.ndarray, control: np.ndarray, penalty: int) -> np.ndarray:
    """
    case = numpy array of read starting positions in the region of interest (case file)
    control = numpy array of read starting positions in the region of interest (control file)
    penalty = integer for factor of lambda to allow for breaking of the inner loop (5 is a good choice)
    Returns a numpy array of floats corresponding to: Start, Stop, Ratio (corresponding to max likelihood)
    
    NOTE: a progress indicator would be great... but cannot be practically implemented in nopython mode.
    """
    lib_size_ratio = 1.0 # Default ratio of both, balances global variation
    lib_size_ratio_calc = (case.size)/(control.size)
    if lib_size_ratio_calc < 0.7 or lib_size_ratio_calc > 1.3:
        # Set ratio lower or higher if it is globally over or under threshold
        lib_size_ratio = lib_size_ratio_calc
    lam = calc_penalty(case.size+ control.size) # Lambda, currently only uses BIC
    M = np.unique(case).size # How many candidate break points are in the case?
    const = penalty * lam # Factor to compare with likelihood to break inner loop
    CBPs = np.unique(case) # Array of candidate break points
    liks = np.zeros(M, dtype=float32) # Initialize array for likelihoods
    # The following calculate the reads on/off for case and control for the initial condition
    case_on = grc_on(CBPs[0], CBPs[0], case)
    case_off = grc_off(CBPs[0], CBPs[0], case)
    control_on = grc_on(CBPs[0], CBPs[0], control)
    control_off = grc_on(CBPs[0], CBPs[0], control)
    # Set initial condition in log likelihood array
    liks[0] = log_lik(case_on, control_on, lib_size_ratio) + \
        log_lik(case_off, control_off, lib_size_ratio) + 0.5 * lam
    first_iteration = True # Used to let the inner loop know when to set the second array item
    results = np.zeros((M-1, 3), dtype=float32)
    for i in range(1, M+1, 1):
        i_pos = CBPs[i] # Position in the block. Minus one to account for 0 based index
        cur_max_l = -9999999999.0 # Set really low, also let numba know it is a float
        cur_max_j = 0 # Start at 0
        res = 0 # Start this over at 0 every loop so that the previous loop doesn't add in
        for j in range(i, 0, -1): 
            # Iterate from beginning of the current break point to the start of the region
            # The goal is to maximize LL and carry forward the best CNV starting point (j_pos)
            j_pos = CBPs[j] 
            case_on = grc_on(j_pos, i_pos, case)
            case_off = grc_off(j_pos, i_pos, case)
            control_on = grc_on(j_pos, i_pos, control)
            control_off = grc_off(j_pos, i_pos, control)
            if first_iteration:
                res = log_lik(case_on, control_on, lib_size_ratio) + \
                    log_lik(case_off, control_off, lib_size_ratio) - lam
                first_iteration = False
            else:
                res = log_lik(case_on, control_on, lib_size_ratio) + \
                    log_lik(case_off, control_off, lib_size_ratio) + liks[j-1] - lam
            if res > cur_max_l: # Is the current likelihood higher than before?
                cur_max_l = res # Set best LL
                cur_max_j = j_pos # Set best start position
            elif res < (cur_max_l - const):
                # If the current res is really low compared to the best, break the loop
                # This provides a significant speed up
                break
        liks[i] = cur_max_l # Set the LL for this i to the best value, this carries (liks[j-1])
        coverage_ratio = ((grc_on(cur_max_j, i_pos, case) + 0.01) / (grc_on(cur_max_j, i_pos, control) + 0.01))
        if coverage_ratio < 1.4 and coverage_ratio > 0.6:
            coverage_ratio = 1.0 
        if coverage_ratio > 5:
            # Eliminate extremes, usually signals that there is no coverage in these areas
            coverage_ratio = 1.0 # Probably a better way to do this
        results[i-1][0] = cur_max_j
        results[i-1][1] = i_pos
        results[i-1][2] = coverage_ratio
    return results

#@jit(float32(int32[:,:], int32[:,:]))
def find_cnv_simple(case, control):
    case_reads = case.size / (np.max(case) - np.min(case))
    control_reads = control.size / (np.max(control) - np.min(control))
    return np.log(case_reads/control_reads)