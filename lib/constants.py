# Natural/Physical Constants

# kT in in kcal/mol*K at 298 K
kT = 0.0019872041 * 298 
# lower boundary for numerical correctness of divisions by zero
num_boundary = 1e-15

def getLogDirPath(WORKDIR, JOBNAME):
    """
    Global definition of the logdir name
    """
    return "{wd}/{jn}-log/".format(wd=WORKDIR, jn=JOBNAME)

def getRunDirPath(WORKDIR, JOBNAME):
    """
    Global definition of the logdir name
    """
    return "{wd}/{jn}-run/".format(wd=WORKDIR, jn=JOBNAME)

