class root_finding:
    
    def n_bisection(self, f, xL, xR, n=100, tolerance=1e-12):
        """
        Uses the bisection method to find values of x between xL and xR for which f(x) = 0, to within the tolerance given.
        Default tolerance is 1e-12, if no tolerance is specified in the function call.
        n is the number of subintervals used to search for roots within the specified interval (default n=100).
        The higher n, the higher computational cost but more likely to find all the roots in the specified interval.
        Returns an array containing all the roots found within the interval.
        """
        # Divide the intervals into n subintervals
        sub_intvls = self.divide_intvl(xL, xR, n) 
        # Initialize the array containg the results of the bisection algorithm for each subinterval
        roots = [] 
        for interval in sub_intvls:
            if f(interval[0])*f(interval[1])<0: 
                roots.append(self.bisection(f, interval[0], interval[1], tolerance)) 
        return roots
    
    
    
    def n_hybrid(self, f, xL, xR, df=None, n=100, tolerance=1e-12):
        """
        Uses a hybrid bisection/newton-raphson method to find values of x between xL and xR for which f(x)=0, to within the tolerance given.
        Default tolerance is 1e-12, if no tolerance is specified in the function call.
        n is the number of subintervals used to search for roots within the specified interval (default n=100).
        The higher n, the higher computational cost but more likely to find all the roots in the specified interval.
        df is the derivative of f. If a derivative function is not specified, the derivative is estimated.
        Returns an array containing all the roots found within the interval
        """
        # Divide the intervals into n subintervals
        sub_intvls = self.__divide_intvl(xL, xR, n)   
        # Initialize the array containg the results of the bisection algorithm for each subinterval
        bisect_ests = []
        for interval in sub_intvls:
            if f(interval[0])*f(interval[1])<0: 
                bisect_ests.append(self.bisection(f, interval[0], interval[1], 0.1))  
                
        # Initialize the array containg the results of the newtown-raphson 
        # algorithm for each bisection starting estimate
        roots = []
        # check whether or not a derivative function was given
        if df==None:
            for est in bisect_ests:
                roots.append(self.newton_raphson(f, est, tolerance=tolerance)) 
        else:
            for est in bisect_ests:
                roots.append(self.newton_raphson(f, est, df, tolerance=tolerance))
            
        return roots
        
            
              
    def bisection(func, xL, xR, tolerance=1e-12):
        """
        Uses the bisection method to find a value of x in (xL,xR) for which f(x)=0, to within the tolerance given.
        Default tolerance is 1e-12 if no tolerance is specified in the function call.
        """
        accuracy = 100 * tolerance
        xM = 0
        
        if func(xL)*func(xR)<0: # make sure the interval has f(x) crossing the y-axis
            while accuracy > tolerance:
                xM = (xL+xR)/2
                if func(xL)*func(xM)<0:
                    xR = xM
                else:
                    xL = xM
                accuracy = abs(func(xM))
            return (xL+xR)/2  # use the midpoint of xL and xR instead of xM (since 
                              # xM is actually equal to either xL or XR at this point)
        else:
            return # the function either has even number of roots in the interval
                   # or has no roots in the interval, in either case return no 
                   # roots found
            
    def newton_raphson(self, f, xinit, df=None, maxiter=1000, tolerance=1e-12):
        """
        Uses the newton-raphson method to find the nearest value of x near xinit for which f(x)=0, to
        within the tolerance given. Default tolerance is 1e-12 if no tolerance is specified in the function
        call.
        """
        itr = 1
        x = xinit
        accuracy = 100*tolerance
        if df==None: #chech whether a derivative function was specified
            while (accuracy > tolerance) & (itr<maxiter):
                x = x - f(x)/self.__derive(f, x) #derivative estimated here
                accuracy = abs(f(x))
                itr = itr + 1
        else:
            while (accuracy > tolerance) & (itr<maxiter):
                x = x - f(x)/df(x)
                accuracy = abs(f(x))
                itr = itr + 1
        return x
        
    ## An estimate of the derivative of an arbitrary 1D function
    def __derive(f, x, h=1e-12):
        """
        Uses the central difference (f(x+h)-f(x-h))/(2 * h) 
        to approximate the value of df/dx where the default
        value for h is h=1e-12
        """
        return (f(x+h) - f(x-h))/(2*h)
        
    
    def __divide_intvl(xL, xR, n):
        """
        Divides the interval (xL, xR) into n subintervals.
        Returns an array of each n subinterval (also stored as an array)\n
        eg. divide(-1,1,2) returns [[-1,0],[0,1]]
        """
        width = (xR - xL)/n
        intvl = []
        for i in range(n):
            intvl.append([xL, xL+width])
            xL = xL+width
        return intvl
