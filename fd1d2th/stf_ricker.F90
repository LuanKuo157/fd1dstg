subroutine stf_ricker(tau,fc,t0,f)
    real:: tau,t0,fc
    real:: u,f0,f,PI=3.141592653
    if (tau<=0.0) then
        f=0.0
    else
       f0=sqrt(PI)/2.0
       u=(tau-t0)*2.0*PI*fc
       f=(u**2.0/4.0-0.5)*exp(-u**2.0/4.0)*f0
   endif
end subroutine
