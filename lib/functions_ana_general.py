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

def weightedHistogram(data, N_bins):
    """
    Returns the normalized histogram of weighted data.
    The data range min(data) to max(data) is divided into N_bins. 
    Takes a data 2d numpy array [[value1, weight1], [value2,weight2], ...] 
    @returns 2d-float [x-values, y-values]
    """
    data = numpy.asarray(data)
    data_min = min(data[:,0])
    data_max = max(data[:,0])
    d         =  float(data_max - data_min ) / N_bins
    hist      =  numpy.zeros([N_bins,2], float)
    # sort data into histogram bins:
    for i in range(0,len(data)):
        index       = int( (data[i,0] - data_min) / d )
        # max(data) entry shall be included in the last bin with index N_bins-1:
        if index==N_bins:
            index = N_bins - 1
        hist[index,1] += data[i,1]
    #Assign the bin positions and normalize:
    total_weight = numpy.sum(data[:,1])
    for i in range(0, N_bins):
        hist[i,0]  = data_min + i * d
        hist[i,1] /= total_weight
   
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

def binIdToCoordinateId(iteration):
    """
    @return An array with indices corresponding between binId and coordinate Id
    """
    sort_indices = []
    for this_bin in iteration:
        sort_indices.append(this_bin.getCoordinateIds())
    return sort_indices

def block_bootstrap(data, function, block_size, number_of_samples = 1000, alpha = 0.05):
    """
    @return Performs a block bootstrap analysis on a time series of sampling data.
    The mean value of a given function of the sampling data and the corresponding confidence intervals are returned.
    """
    from random import randint
    
    # generate a resampled dataset based on non-overlapping blocks
    def resample(data, block_size):
        number_of_blocks = int (len(data) / block_size )
        data_resampled = []
        for j in range(0,number_of_blocks):
            random_block_index = randint(0, number_of_blocks - 1)
            data_block_tmp = data[block_size * random_block_index:block_size * (random_block_index + 1) ]
            data_resampled = numpy.append(data_resampled, data_block_tmp)
  
        return data_resampled

    # calculate the alpha percentiles for lower and upper boundaries
    def percentile(data, alpha):
        data     = numpy.sort(data)
        for i in range(len(data)):
            if i >= len(data) * alpha:
                p_l = data[i]
                break       
        for i in range(len(data)-1 ,0-1,-1):
            if i >= len(data) * alpha:
                p_u = data[i]      
                break
        return p_l, p_u
            
    # Evaluate function on the resampled datasets
    function_values = []
    for i in range(0, number_of_samples):
        data_resampled = resample(data, block_size)
        function_values.append( function(data_resampled) )
    
    return ( numpy.mean(data), percentile(function_values, alpha) )
    
        
    print numpy.mean(function_values), numpy.std(function_values)
    
def cumulative_mean(data):
    """
    @return cumulative mean of data
    """
    cumulative_mean = []
    cumulative_sum = 0.0        
    for i, value in enumerate(data):
        cumulative_sum += value
        cumulative_mean.append(cumulative_sum/float(i+1))
    return cumulative_mean
        
def getStateFromCoordinate(segment, state_A, state_B):
    """
    Returns the state of a segment
    @return string state
    """
    state_per_dimension = []
    # lazy 1d implementation
    for coordinate in [segment.getCoordinates()[0]]:
        if coordinate > state_A[0] and coordinate < state_A[1]:
            state_per_dimension.append('A')
        elif coordinate > state_B[0] and coordinate < state_B[1]:
            state_per_dimension.append('B')
        else:
            state_per_dimension.append('0')
    
    this_state = state_per_dimension[0]
    for state in state_per_dimension[1:]:
            if state != this_state:
                this_state = '0'
                break 
    return this_state   
    
    
    