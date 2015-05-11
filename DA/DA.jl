module DA
export threeDVar,EKF,EnKF,ETKF

EKF = function(xf::Array,y::Array,H::Array,R::Array,pf::Array)
    # compute the Kalman gain matrix
    K = \( *(pf,H.') , R+ *(H ,*(pf,H.') ) );

    # innovation
    d = y - *(H,xf);

    # analysis
    xa = xf + *(K,d);
    pa = *( (eye(length(xf)) - *(K,H)), pf );

    return xa,pa
end

threeDVar = function(xf::Array,y::Array,H::Array,R::Array,B::Array)
    # compute the intermediary matrix solution
    K = \( ( inv(B) + *(H.', *(inv(R),H)), *(H.',inv(R) ) ) );

    # innovation
    d = y - *(H,xf);

    # analysis
    xa = xf + *(K,d); 

    return xa
end

EnKF = function(X_f::Array,y_0::Array,H::Array,R::Array,delta::Float64)
    # gerneralized EnKF
    # "perturbed observation" method
    #
    # INPUTS
    # X_f: array of the forecast states
    #      forecasts going down, models going across
    #      X_f[:,1] should be the first forecast
    # y_0: observation
    # H: observation operator
    # R: observation error covariance
    # delta: multiplicative inflation 
    #        for analysis ensemble

    N = size(X_f,2)

    # set the error the randomly perturbing the analysis IC for ensemble members
    # this is the std dev of a normal error
    pertErr = R[1,1]

    # let x_f now be the average
    x_f = mean(X_f,2)

    # this computes the difference from the mean without that forecast
    X_f_diff = X_f - repmat(x_f,1,N)
    println("X_f_diff:")
    println(size(X_f_diff))
    println(X_f_diff[1:10,:])

    # multiplicative inflation
    X_f_diff = sqrt(1+delta).*X_f_diff

    # estimate error covariance
    pf = 1/(N-1).*( *(X_f_diff,X_f_diff.') )
    println(pf[1:10,1:10])
    
    # compute Kalman gain matrix
    K = \( *(pf,H.') , R+ *(H ,*(pf,H.') ) )
    println(K[1:10,1:10])

    # innovation
    d = y_0 - *(H,x_f)
    println(d[1:10])

    X_a = zeros(size(X_f))
    for j=1:N
        pertY = y_0+randn(size(y_0)).*pertErr
    
        # innovation
        d = pertY - *(H,(x_f+X_f_diff[:,j]))

        X_a[:,j] = x_f+X_f_diff[:,j] + K*d
    end

    return X_a
end

ETKF = function(X_f::Array,y_0::Array,H::Array,R::Array,delta::Float64)
    # Ensemble Transform Kalman Filter
    #
    # INPUTS
    # X_f: array of the forecast states
    #      forecasts going down, models going across
    #      X_f[:,1] should be the first forecast
    # y_0: observation
    # H: observation operator
    # R: observation error covariance
    # delta: multiplicative inflation 
    #        for analysis ensemble

    N = size(X_f,2)

    # set the error the randomly perturbing the analysis IC for ensemble members
    # this is the std dev of a normal error
    pertErr = R[1,1]

    # let x_f now be the average
    x_f = mean(X_f,2)

    # this computes the difference from the mean without that forecast
    X_f_diff = X_f - repmat(x_f,1,N)
    println(X_f_diff[1:10,:])

    # multiplicative inflation
    X_f_diff = sqrt(1+delta).*X_f_diff

    # estimate error covariance
    pf = 1/(N-1).*( *(X_f_diff,X_f_diff.') )
    println(pf[1:10,1:10])
    
    # compute Kalman gain matrix
    K = \( *(pf,H.') , R+ *(H ,*(pf,H.') ) )
    println(K[1:10,1:10])

    # innovation
    d = y_0 - *(H,x_f)
    println(d[1:10])

    x_a = x_f + *(K,d)

    # perform the transform 
    # it is (from the literature)
    # absolutely necessary to perform this inversion
    # D = *(H,X_f_diff)
    # C = *(R,*(H,X_f_diff))
    # B = \( D.' , C )
    # A = (N-1).*eye(N)+B
    # p_a = inv(A)
    # p_a = inv( (N-1).*eye(N)+ *( \( *(H,X_f_diff).' , R) ,*(H,X_f_diff) ) )
    # from matlab:
    # p_a = inv((N-1)*eye(N)+(H*X_f_diff)'/(R)*(H*X_f_diff));
    C = (N-1).*eye(N)
    B = *( H , X_f_diff )
    # println(size(B))
    # println(size(R))
    # this bugger is tricky
    # julia polyalgorithm reshapes the output, so don't transpose input
    D = \( B , R )
    A = C+*( D , B)
    p_a = inv(A)
    # println("analysis error covariance:")
    # println(typeof(p_a))
    # println(issym(p_a))
    # println("difference from symmetry of first row:")
    # println(p_a[:,1]-p_a[1,:].')
    # (the difference is on the order of machine epsilon, for simple test)
    # so force the symmetry

    T = sqrtm(Symmetric((N-1).*p_a))
    # println("transform type is:")
    # println(typeof(T))

    X_a_diff = *(X_f_diff,T)

    # now compute the analysis for each ensemble member
    X_a=repmat(x_a,1,N) + X_a_diff

    return X_a
end

end
