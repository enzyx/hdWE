import numpy

def histogram(data, N_bins):
    """
    Returns the normalized histogram of data.
    The data range min(data) to max(data) is divided into N_bins. 
    """
    data_min = min(data)
    data_max = max(data)
    d         =  float(data_max - data_min ) / N_bins
    hist      =  numpy.zeros([N_bins,2], float)
    # sort data into histogram bins:
    for i in range(0,len(data)):
        index       = int( (data[i] - data_min) / d )
        # max(data) entry shall be included in the last bin with index N_bins-1:
        if index==N_bins:
            index = N_bins - 1
        hist[index,1] += 1
    #Assign the bin positions and normalize:
    for i in range(0, N_bins):
        hist[i,0]  = data_min + i * d
        hist[i,1] /= len(hist)
   
    return hist
    
def autocorrelation_function(x):
    """
    Returns the autocorrelation function of a data set, using numpy.corrcoef.
    """
    autocorrelation_function = []
    for t in range(len(x)/2):
        correlation_coefficient = numpy.corrcoef(numpy.array([x[0:len(x)-t], x[t:len(x)]]))[0,1] 
        autocorrelation_function.append(correlation_coefficient)
    return autocorrelation_function 