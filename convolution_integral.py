import pandas as pd
import numpy
import scipy
import pylab
import pdb
from scipy.interpolate import interp1d
from scipy.integrate import trapz
from matplotlib.dates import date2num
import datetime
import lmfit as lm

def find(condition):
    res, = np.nonzero(np.ravel(condition))
    return res

def make_conv_ts_tau(C_t,tracer,t_obs,tau_vec,lamba,disp_fac=.5,mod_type='exponential',trit_flag=0):
    '''Returns the vector of possible conentrations of tracer at observation time t_obs for a given age mixing type (mod_type), at
    a given time of observation (t_obs) for a range of mean ages (tau vec), given the historical input data C_t from a time index pandas dataframe'''
    C_conv_tau = numpy.zeros(len(tau_vec));
    for i in range(len(tau_vec)):
        tau_i = tau_vec[i];
        C_tau_i = conv_int_discrete(C_t,tracer,t_obs,tau_i,lamba,disp_fac=disp_fac,mod_type=mod_type,trit_flag=trit_flag);
        C_conv_tau[i] = C_tau_i;
    return C_conv_tau

def make_conv_ts(C_t,tracer,t_vec,tau,lamba,disp_fac=.5,mod_type='exponential',trit_flag=0):
    """ creates a time series of the expected concentration at each time in t_vec for a given aquifer and mixing model, with given input time series
    C_t - the conentration for each time in t_vec. Returns C_conv, the convoloved concentration at each time point in t_vec."""
    C_conv = numpy.zeros(len(t_vec));
    print('Calculating the convolved time series for tau =' + str(tau/365) + ' years')
    for i in range(len(t_vec)):
        t_i = t_vec[i];
        C_t_i = conv_int_discrete(C_t,tracer,t_i,tau,lamba,disp_fac=disp_fac,mod_type=mod_type,trit_flag=trit_flag);
        C_conv[i]=C_t_i;
    return C_conv;
    
def conv_int_discrete(C_t,tracer,t_i,tau,lamba,disp_fac=.5,mod_type='exponential',trit_flag=0):
    """ Convulution integral program calculates the the concentration at observation time t_i (numpy datetime of the observation time)  as the convolution integral for a given model type. Needs the input concentration time series, C_t is a pandas data frame indexed by time with historical input, for example GNIP.data, tracer is the key for the tracer to be convolved, t_i (the observation time ), and aquifer mean age.  Returns the concentration at time t_i for the mean residence time tau."""
    #extend the data serie back in time
    #pdb.set_trace()
    C_t_start = C_t[:C_t.index[0]]
    padding = int(tau*10) #for now will keep daily resolution, but this could an issue for long residence times.
    ix = pd.Index(numpy.arange(padding))
    C_t_start = C_t_start.reindex(ix)
    if tracer in ('O18','H2'):
        C_t_start.iloc[-1]=C_t.mean().values
    else:
        C_t_start.iloc[-1]=C_t.iloc[0].values
        
    C_t_start = C_t_start.fillna(method='bfill')
    #only use concentration history before the sampling date
    #C_t = C_t[:t_i]
    if C_t.index.freqstr != 'D':
        C_t_resampled = C_t.resample('D').interpolate(method='time',inplace=False) #make sure we're at daily resolution
    else:
        C_t_resampled = C_t
    
    #calculate the length of the input record before the observation time in days    
    obs_lag = t_i-C_t.index[0]  

    #now only use concentration history before the sampling date
    c_i = C_t_resampled[:t_i][tracer]
    
    #calculate C(t-tau)
    lag_times = (t_i-c_i.index).days.values/1. #calculate the days (need to be a float)
    lag_times[lag_times[:]==0]=.0001
    #lag_times[0]=.0001
    #delta_t = t_vec_i[1]-t_vec_i[0];
    
    if mod_type == 'piston':
        g_t = numpy.zeros(len(lag_times));
        ix = find(abs(lag_times-tau)==min(abs(lag_times-tau)));
        g_t[ix]=1.;

    if mod_type == 'exponential':  #uniform recharge on an unconfined aquifer.
        g_t = (1./tau)*numpy.exp(-lag_times/tau);
    
    if mod_type == 'linear':  #uniform recharge on a linearly increasing with depth aquifer
        g_t = numpy.zeros(len(lag_times));
        for i in range(len(lag_times)):
            tt_i = lag_times[i];
            if tt_i <= 2*tau:
                g_t[i] = 1./(2*tau);
                
    if mod_type == 'exp_pist_flow': #recharge area x<x_star, confined of constant thickness after that
        #x_star = 100.; #2km flow path length
        #x = 100.; #recharge area
        eta = 1.5; #eta is volume total/volume of exponential aquifer
        g_t = numpy.zeros(len(lag_times));
        #tau_mean = (tau*(x+x_star))/(x);
        for i in range(len(lag_times)):
            #if i == 0:
             #   if error_flag = 0:
             #       print 'warning volume fraction eta hard coded at this point'
            #      error_flag = 1;
            tt_i = lag_times[i];
            if tt_i > tau*(1-(1/eta)):
                g_t[i] = (eta/tau)*numpy.exp(-eta*tt_i/tau+eta-1.);
    
    if mod_type == 'lin_pist_flow': #recharge x<x_star, constant thickness in confined part
        #x_star = 100; #2km flow path length
        #x = 100; #recharge area
        g_t = numpy.zeros(len(lag_times));
        #tau_mean = (tau*(x+x_star))/(x);
        eta = 1.5        
        for i in range(len(lag_times)):
           # if i == [0]:
            #    if error_flag = 0:
             #       print 'warning volume fraction eta hard coded at this point';
             #       error_flag = 1;
            tt_i = lag_times[i];
            if tt_i >= tau - tau/eta:
                if tt_i <= tau + tau/eta:
                    g_t[i] = eta/(2.*tau);
    
    if mod_type == 'dispersion':
        g_t = ((1./tau)/(numpy.sqrt(4.*numpy.pi*disp_fac*(lag_times/tau))))*(1./(lag_times/tau))*numpy.exp(-1.*(((1.-(lag_times/tau))**2)/(4.*disp_fac*(lag_times/tau))));
        # disp = dispersion parameter = 1/Pe = D/vx
    
    #some specific half lives...
    if tracer == 'H3':
        t_half = 4500 #days from Lucas 2000
        lamba = -1*numpy.log(0.5)/t_half

    #add in exponential decay
    #G_t = g_t*numpy.exp(-lamba*t_prime)

    #accumulate if tritiogenic helium
    if tracer == 'He3t':
        t_half = 4500 #days from Lucas 2000
        lamba = -1*numpy.log(0.5)/t_half
    #    G_t = g_t*(1-numpy.exp(-lamba*t_prime))
    
    # convolve - I think this is where a massive speed up could occur.
    C_t_i = 0;
    for k in range(len(lag_times)):
        if mod_type=='piston':
            C_step = c_i[k]*g_t[k]*numpy.exp(-lamba*lag_times[k]);
        else:    
            if lag_times[k] == .0001:
                lag_time1 = lag_times[k-1]/2.;
                lag_times[k] = lag_time1;
                C_step =c_i[k]*(g_t[k])*numpy.exp(-lamba*lag_times[k])*lag_time1
        #if mod_type == 'piston':
        #    C_step = c_i[k]*g_t[k]*numpy.exp(-lamba*lag_times[k]);
            else:
                try:
                    C_step =c_i[k]*g_t[k]*numpy.exp(-lamba*lag_times[k])*(lag_times[k]-lag_times[k+1]);
                except IndexError:
                    C_step =c_i[k]*g_t[k]*numpy.exp(-lamba*lag_times[k])*(lag_times[k]-0);
        C_t_i = C_t_i + C_step;
        
        
    print(t_i,C_t_i)
    return C_t_i;


def conv_int_fft(C_t,tracer,t_i,tau,lamba,disp_fac=.5,mod_type='exponential',trit_flag=0,bbar=0.0001,Phi_im=0.01):
    """ Convolution integral program calculates the the concentration at observation time t_i (numpy datetime of the observation time)  as the convolution integral for a given model type. Needs the input concentration time series, C_t is a pandas data frame indexed by time with historical input, for example GNIP.data, t_i (the observation time ), and aquifer mean age.  Returns the concentration history at the site ending at time t_i for the mean residence time tau.  Uses numpy.convolve and is WAY faster than the discrete method.  Only use discrete as a final checking method."""

    #extend the data serie back in time
    C_t_start = C_t[:C_t.index[0]]
    padding = int(tau*10) #for now will keep daily resolution, but this could an issue for long residence times.
    ix = pd.Index(numpy.arange(padding))
    C_t_start = C_t_start.reindex(ix)
    if tracer in ('O18','H2'):
        C_t_start.iloc[-1]=C_t.mean().values
    else:
        C_t_start.iloc[-1]=C_t.iloc[0].values
        
    C_t_start = C_t_start.fillna(method='bfill')
    #only use concentration history before the sampling date
    #C_t = C_t[:t_i]
    #C_t_resampled = C_t.resample('D').sum().interpolate() #make sure we're at daily resolution everywhere - the sum is somehow needed and deals with time offset I think
    #C_t_resampled = C_t_resampled.interpolate() # linear interpolation
    if C_t.index.freqstr != 'D':
        C_t_resampled = C_t.resample('D').interpolate(method='time',inplace=False) #make sure we're at daily resolution
    else:
        C_t_resampled = C_t
    ix2 = numpy.arange(len(C_t_start),len(C_t_start)+len(C_t_resampled)) #new int day index
    C_t_resampled['iix'] = ix2
    C_t_int = C_t_resampled.set_index(ix2)
    del C_t_int['iix']
    C_t_extended = pd.concat([C_t_start,C_t_int],axis=0)
    c_i = C_t_extended[tracer]
    
    #g(t')
    t_prime = numpy.arange(0,10*tau) #times for calculating the weighting function
    t_prime[t_prime[:]==0]=.0001
    #lag_times[0]=.0001
    #delta_t = t_vec_i[1]-t_vec_i[0];
    
    if mod_type == 'piston':
        g_t = numpy.zeros(len(t_prime));
        ix = find(abs(t_prime-tau)==min(abs(t_prime-tau)));
        g_t[ix]=1.;

    if mod_type == 'exponential':  #uniform recharge on an unconfined aquifer.
        g_t = (1./tau)*numpy.exp(-t_prime/tau);
    
    if mod_type == 'linear':  #uniform recharge on a linearly increasing with depth aquifer
        g_t = numpy.zeros(len(t_prime));
        for i in range(len(t_prime)):
            tt_i = t_prime[i];
            if tt_i <= 2*tau:
                g_t[i] = 1./(2*tau);
                
    if mod_type == 'exp_pist_flow': #recharge area x<x_star, confined of constant thickness after that
        #x_star = 100.; #2km flow path length
        #x = 100.; #recharge area
        eta = 1.5; #eta is volume total/volume of exponential aquifer
        g_t = numpy.zeros(len(t_prime));
        #tau_mean = (tau*(x+x_star))/(x);
        for i in range(len(t_prime)):
            #if i == 0:
             #   if error_flag = 0:
             #       print 'warning volume fraction eta hard coded at this point'
            #      error_flag = 1;
            tt_i = t_prime[i];
            if tt_i > tau*(1-(1/eta)):
                g_t[i] = (eta/tau)*numpy.exp(-eta*tt_i/tau+eta-1.);
    
    if mod_type == 'lin_pist_flow': #recharge x<x_star, constant thickness in confined part
        #x_star = 100; #2km flow path length
        #x = 100; #recharge area
        g_t = numpy.zeros(len(t_prime));
        #tau_mean = (tau*(x+x_star))/(x);
        eta = 1.5        
        for i in range(len(t_prime)):
           # if i == [0]:
            #    if error_flag = 0:
             #       print 'warning volume fraction eta hard coded at this point';
             #       error_flag = 1;
            tt_i = t_prime[i];
            if tt_i >= tau - tau/eta:
                if tt_i <= tau + tau/eta:
                    g_t[i] = eta/(2.*tau);
    
    if mod_type == 'dispersion':
        g_t = ((1./tau)/(numpy.sqrt(4.*numpy.pi*disp_fac*(t_prime/tau))))*(1./(t_prime/tau))*numpy.exp(-1.*(((1.-(t_prime/tau))**2)/(4.*disp_fac*(t_prime/tau))));
        # disp = dispersion parameter = 1/Pe = D/vx
        
    if mod_type == 'frac_inf_diff':
        #modified from painter 2008 and villermeax 1981
        D_o = (2.3e-9)*60*60*24 #diffusion coefficient in water m2/s (time is all in days...)
        #Phi_im = 0.01 #immobile zone porosity
        R_im = 1 #immoble zone retardation for linear sorption
        kappa = Phi_im*numpy.sqrt((D_o*Phi_im**2)*R_im)
        tt_prime = t_prime # the times to solve for the diffusion modified transit distribution
        f_t_tran = numpy.zeros(len(tt_prime))
        Beta_bar = tau/bbar
        for i in range(len(tt_prime)):
            t_tran = tt_prime[i]
            # calculate the advective time distribution up to the transit time distribution
            tadv = numpy.logspace(-6,numpy.log10(t_tran-1e-6),num=1000) # this discretization might be overkill... #trying to get rid of the singularity...
            f_tadv = ((1./tau)/(numpy.sqrt(4.*numpy.pi*disp_fac*(tadv/tau))))*(1./(tadv/tau))*numpy.exp(-1.*(((1.-(tadv/tau))**2)/(4.*disp_fac*(tadv/tau)))) # first the advective travel time disribution  - dispersive only for now
            
            t_ret = t_tran - tadv
            Beta_tadv = tadv/bbar
            f_ret_Beta_tadv = (kappa*Beta_tadv)/(2*numpy.sqrt(numpy.pi)*t_ret**(3/2))*numpy.exp((-1*kappa**2*Beta_tadv**2)/(4*t_ret)) # assumes uniform fracture aperature over the transport length.
            f_ret_Beta_tadv = f_ret_Beta_tadv/trapz(f_ret_Beta_tadv[::-1],t_ret[::-1]) #normalize - this is key and was a big stumbling block - should it be published?
            f_i = f_ret_Beta_tadv*f_tadv
            f_t_tran[i] = trapz(f_i,tadv)
        
        # for interpolation reasons...
        f_t_tran[0]=0
        tt_prime[0]=0
        res2 = trapz(f_t_tran,tt_prime)
        f_t_tran = f_t_tran/res2 #normalize distribution.
        #resample at t_prime:
        f2 = interp1d(tt_prime,f_t_tran,kind="linear")
        g_t = f2(t_prime)
        #pdb.set_trace()
        mean_travel_time = trapz(g_t*t_prime,t_prime)/trapz(g_t,t_prime)/365.
        print("mean travel time is " + '%3.2f' %mean_travel_time)

    #some specific half lives...
    if tracer == 'H3':
        t_half = 4500 #days from Lucas 2000
        lamba = -1*numpy.log(0.5)/t_half

    #add in exponential decay
    G_t = g_t*numpy.exp(-lamba*t_prime)

    #accumulate if tritiogenic helium
    if tracer == 'He3t':
        t_half = 4500 #days from Lucas 2000
        lamba = -1*numpy.log(0.5)/t_half
        G_t = g_t*(1-numpy.exp(-lamba*t_prime))
    
    C_t_c = numpy.convolve(c_i,G_t,mode='valid')
    dt = pd.Timedelta(len(C_t_c)-1,'D')
    tstart = t_i - dt
    ix = pd.DatetimeIndex(start=tstart,end=t_i,freq='D')
    C_t_c = pd.Series(C_t_c,index=ix)
    print('mean age = %1.2g' %(tau/365.))
    print('Conctration at ', t_i, '= %1.3g' %C_t_c[-1])
        
    return C_t_c

def residual(params,C_t,tracer,errror_perc,t_i,df_gw,lamba=0.0,mod_type='exponential',show_interim_results=True):
    '''calulates the least squares terms for a convolved times series and a groundwater data set for the given tracer, returns the array of least square terms.'''
    tau = params['mean_age'].value
    disp_fac = params['disp_fac'].value
    c_conv = conv_int_fft(C_t,tracer,t_i,tau,lamba,disp_fac=disp_fac,mod_type=mod_type,trit_flag=0.)
    df_gw_res = df_gw[:c_conv.index[-1]][tracer]
    df_gw_res = df_gw_res.dropna()
    c_conv_rs = c_conv.reindex(df_gw_res.index)
    residual = (df_gw_res-c_conv_rs)/(c_conv_rs*(error_perc/100.))
    if show_interim_results:
        print('tau = ' + str(tau/365) + ' , D = '+ str(disp_fac), 'res = ' + str(numpy.linalg.norm(residual)))
        pylab.figure()
        pylab.plot(c_conv_rs.index,c_conv_rs,'b-',label=r'$\tau =$' + str(tau/365) + r'  $ \prime{\bar{D}}$ = ' + str(disp_fac))
        pylab.plot(df_gw_res.index,df_gw_res,'ro',label='data')
        #pylab.plot(residual.index,residual,'g-',label='residual')
        pylab.legend(loc='best')
        pylab.show()
    return residual.values

def residualmt(params,C_t,tracers,error_perc,t_i,df_gw,lamba=0.0,mod_type='exponential',show_interim_results=True):
    '''calulates the least squares terms for a convolved times series and a groundwater data set for tracers in tracer list, returns the array of least square terms.'''
    tau = params['mean_age'].value
    disp_fac = params['disp_fac'].value
    residual = numpy.array([])
    for tracer in tracers:
        c_conv = conv_int_fft(C_t,tracer,t_i,tau,lamba,disp_fac=disp_fac,mod_type=mod_type,trit_flag=0.)
        #counting on c_conv being daily here, and daily sampled.  This will probably break some day.
        df_gw_res = df_gw[:c_conv.index[-1]][tracer]
        df_gw_res = df_gw_res.dropna()
        c_conv_rs = c_conv.reindex(df_gw_res.index)
        residual_i = (df_gw_res-c_conv_rs)/(c_conv_rs*(error_perc[tracer]/100.))
        residual = numpy.concatenate([residual,residual_i])
        if show_interim_results:
            print('tau = ' + str(tau/365) + ' , D = '+ str(disp_fac), 'res = ' + str(numpy.linalg.norm(residual)))
            pylab.figure()
            pylab.plot(c_conv_rs.index,c_conv_rs,'b-',label=r'$\tau =$' + str(tau/365) + r'  $ \prime{\bar{D}}$ = ' + str(disp_fac))
            pylab.plot(df_gw_res.index,df_gw_res,'ro',label='data')
            #pylab.plot(residual.index,residual,'g-',label='residual')
            pylab.legend(loc='best')
            pylab.show()
    return residual

def estimate_mean_age(tau_i,tau_high,tau_low,disp_fac,disp_fac_high,disp_fac_low,C_t,tracer,error_perc,t_i,df_gw,lamba=0.0,mod_type='exponential',est_disp=False,show_interim_results=True):
    '''runs lmfit for a non-linear least squares fit of lumped parameter model to observed data. needs the initial mean age guess, bounds on the mean age, the initial dispersion factor (only used for dispersion model, ignored otherwise), and bounds on the dispersion factor, the inpute concentration time series as a pandas time indexed data frame, the tracer being matched (name must be the same in input and observed data), an estimate of the percent error of the data, the observed data as a pandas time indexed data frame, the decay coefficient for the tracer, the lumped parameter model type, and a boolean indicating weather or not to solve for the dispersion factor (dispersion model only!).  Returns the lm.minimizerresult obj with all the data on the fit.  Uses Levenberg-Marquahrt non-linear least squares for the fit. '''
    x = lm.Parameters()
    x.add('mean_age',value=tau_i,min=tau_low,max=tau_high,vary=True)
    if mod_type != 'dispersion':
        est_disp = False
    x.add('disp_fac',value=disp_fac,min=disp_fac_low,max=disp_fac_high,vary=est_disp)
    #mini = lm.Minimizer(residual,x,fcn_args=(C_t,tracer,t_i,df_gw),fcn_kws={'lamba':lamba,'disp_fac':disp_fac,'mod_type':mod_type},scale_covar=False,ftol=1e-5,xtol=1e-5,epsfcn=.001)
    #out = mini.minimize()
    #pdb.set_trace()
    if type(tracer) == str:
        out = lm.minimize(residual,x,args=(C_t,tracer,error_perc,t_i,df_gw),kws={'lamba':lamba,'mod_type':mod_type,'show_interim_results':show_interim_results},scale_covar=False,ftol=1e-5,xtol=1e-5,epsfcn=.001)
    else:
        out = lm.minimize(residualmt,x,args=(C_t,tracer,error_perc,t_i,df_gw),kws={'lamba':lamba,'mod_type':mod_type,'show_interim_results':show_interim_results},scale_covar=False,ftol=1e-5,xtol=1e-5,epsfcn=.001)
    print((lm.fit_report(out)))
    return out
