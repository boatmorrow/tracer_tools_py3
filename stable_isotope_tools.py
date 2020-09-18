import pdb
import pandas as pd
import numpy as N
import datetime

#module containing functionality for working with stable and radio isotopes of water
class GNIP:
    '''GNIP class for working with downloaded GNIP data'''
    def __init__(self,name):
        self.name=name
        self.data=pd.DataFrame()
        self.sitename=''
        self.latitude=0.
        self.longitude=0.
        self.elevation=0.
        self.alpha=0.

    #read a GNIP spreadsheet
    def read_gnip(self,ifile):
        '''reads in a GNIP excel spreadsheed and fills out data for the GNIP class.  Timeseries data is a pandas dataframe'''
        df = pd.read_excel(ifile,usecols="M,Q,S,U,X,Y")
        df['Date'] = pd.to_datetime(df['Date'])
        #dt = pd.Timedelta('15 days')
        #pdb.set_trace()
        #for i in xrange(len(df)):
        #    df['Date'].ix[i]=df['Date'].ix[i]+dt
        #df = df.set_index('Date')
        #self.data=df
        #get site info
        df2 = pd.read_excel(ifile,header=0,parse_cols="A:F",index_col=None)
        self.sitename = df2['Site'][0]
        self.latitude = df2['Latitude'][0]
        self.longitude = df2['Longitude'][0]
        self.elevation = df2['Altitude'][0]

#    def calc_alpha(self,gw_mean_conc,tracer='H2'):
#        '''calculate alpha = alpha_s/alpla_w where alpha_i=I/P is the runoff ratio, using equation 2.15 from the IAEA book.  Needs the mean gw concentration measured in an area near the GNIP site. Tracer should be either 'H2' or 'O18'. For now assumes Northern hemisphere, 2 6 month seasons. '''
#        #process the precip
#        #make a seasonal period_range - assumes summer starts in march
#        df_b = datetime.datetime(self.data.index[0].year-1,3,1)
#        df_e = self.data.index[-1]
#        #add the precip weighted isotopes
#        self.data.loc[:,'pw_tracer']=self.data[tracer]*self.data['Precipitation'].copy()
#        prng = pd.period_range(df_b,df_e,freq='6M')
#        #keep only full seasons for calculating alpha = for now only two six month seasons, this could be vastly improved
#        df = self.data.groupby(prng.asof).filter( lambda x: len(x) == 6)
#        df_grouped = df.groupby(prng.asof)
#        df_seasonal_sum = df_grouped.sum()
#        #now inner join these two to process the the two
#        df_ac = df_seasonal_sum[['Precipitation','pw_tracer']]
#        #keep only full years of data
#        prng = pd.period_range(datetime.datetime(df_ac.index[0].year-1,3,1),datetime.datetime(df_ac.index[-1].year+1,3,1),freq='12 M')
#        df_ac = df_ac.groupby(prng.asof).filter(lambda x: len(x) == 2)
#        #equation 2.15 IAEA isotope book
#        ict = 0
#        wflag = 0
#        sflag = 0
#        last_year = 0
#        SpdT_w = 0
#        dSp_w = 0
#        SpdT_s = 0
#        dSp_s = 0
#        dT_gw = gw_mean_conc
#        for y in df_ac.index[:].year:
#            if y != last_year:
#                wflag = 0
#                sflag = 0
#                last_year = y
#            if df_ac.index[ict].month == 9:
#                wflag = 1
#                SpdT_w = SpdT_w + df_ac['pw_tracer'][ict]
#                dSp_w = dSp_w + dT_gw*df_ac['Precipitation'][ict]
#            if df_ac.index[ict].month == 3:
#                sflag = 1
#                SpdT_s = SpdT_s + df_ac['pw_tracer'][ict]
#                dSp_s = dSp_s + dT_gw*df_ac['Precipitation'][ict]
#                last_year = y
#            ict += 1
#        alpha = (SpdT_w - dSp_w)/(dSp_s-SpdT_s)
#        self.alpha = alpha
#        return alpha
    def read_gnip_csv(self,ifile):
        '''reads in a GNIP csv and fills out data for the GNIP class.  Timeseries data is a pandas dataframe'''
        df = pd.read_csv(ifile,usecols=[12,16,18,20,23,24])
        df['Date'] = pd.to_datetime(df['Date'])
        #dt = pd.Timedelta('15 days')
        #pdb.set_trace()
        #for i in xrange(len(df)):
        #    df['Date'].ix[i]=df['Date'].ix[i]+dt
        df = df.set_index('Date')
        self.data=df
        #get site info
        df2 = pd.read_csv(ifile,header=0,usecols=[2,5,6,7],index_col=None)
        self.sitename = df2['Site'][0]
        self.latitude = df2['Latitude'][0]
        self.longitude = df2['Longitude'][0]
        self.elevation = df2['Altitude'][0]

def weighted_iso_ts(precip_ts,gw_ts,alpha,tracer,iw_f='ps_func'):
    '''make an infiltration weighted time series according to equation 2.18 IAEA book ('iaea') or my modified and better equation 'ps_func', needs a preciptation pandas time index dataframe, a groundwater time indexed data frame with the the same tracer names, and the alpha=alpha_s/alpha_w partition coefficient.  Assumes northern hemisphere and two 6 month seasons. Return a data frame with all the precip data with an Iw_tracer column added which is the infiltration weighted infiltration isotope value.  Annual infil weighting only for tritium.  Suggest using an average alpha for H2 and O18.'''
    print('average precip H2 = ' + str(precip_ts['H2'].mean()))
    print('average gw H2 = ' + str(gw_ts['H2'].mean()))
    print('average precip O18= ' + str(precip_ts['O18'].mean()))
    print('average tunnel O18= ' + str(gw_ts['O18'].mean()))
    
    # if we're doing this for tritium - annunual infil weighting only.
    if tracer == 'H3':
        iw_f = 'an_av'

    #calculate yearly infiltration
    #first keep only full years of data - with seasonal starts
    prng = pd.period_range(datetime.datetime(precip_ts.index[0].year-1,3,1),datetime.datetime(precip_ts.index[-1].year+1,3,1),freq='12 M')
    df = precip_ts.groupby(prng.asof).filter(lambda x: len(x) == 12)
    #add alpha_i where alpha_w = 1 and alpha_s = alpha
    df.loc[:,'alpha']=N.ones(len(df.index))
    for k in df.index:
        if k.month >= 3 and k.month <9:
            df.loc[k]['alpha']=alpha
    #monthly infiltration
    df.loc[:,'infil']=df['alpha']*df['Precipitation']
    df.loc[:,'infil_tracer']=df['infil']*df[tracer]

    #now group by seasons and calulate annual totals
    prng = pd.date_range(datetime.datetime(df.index[0].year-1,3,1),datetime.datetime(df.index[-1].year+1,3,1),freq='12M')
    df_seasonal = df.groupby(prng.asof)
    df_seasonal = df_seasonal.sum()
    df_an_av = df_seasonal['infil_tracer']/df_seasonal['infil']
    dtt = datetime.datetime(df_an_av.index[0].year,df_an_av.index[0].month,df_an_av.index[0].day)-datetime.datetime(df_an_av.index[0].year,df_an_av.index[0].month,df.index[0].day) #make days match
    df_an_av = pd.DataFrame(df_an_av)
    df_an_av.columns=['an_av_iw']
    ix = df_an_av.index-dtt
    df_an_av = df_an_av.set_index(ix)
    df_an_av = df_an_av.reindex(df.index)
    df_an_av = df_an_av.interpolate()
    #pdb.set_trace()
    
    # equation 2.18 (IAEA book)
    df.loc[:,'Iw_'+tracer]=N.zeros(len(df.index))
    bar_d = df[tracer].mean()  #not sure weather bar_d should be the average of precip or groundwater... But the more I think about it, it should be groundwater...
    bar_d = gw_ts[tracer].mean()
    if iw_f == 'an_av':
        df['Iw_'+tracer]=df_an_av['an_av_iw']
    else:
        for k in df.index:
            try:
                if iw_f == 'iaea':
                    df.loc[k]['Iw_'+tracer] = bar_d + (df.loc[k][tracer] - bar_d)* df.loc[k]['infil']/(df_seasonal[df_seasonal.index.year==k.year].infil/12)
                else:
                    df.loc[k]['Iw_'+tracer] = bar_d + (df.loc[k][tracer] - bar_d)* df.loc[k]['infil']/(df_seasonal[df_seasonal.index.year==k.year].infil) #no longer 2.18, but something I think makes more sense.
            except ValueError:
                if k.year > df_seasonal.index[-1].year and k.month < 3:
                    if iw_f == 'iaea':
                        df.loc[k]['Iw_'+tracer] = bar_d + (df.loc[k][tracer] - bar_d)* df.loc[k]['infil']/(df_seasonal.loc[df_seasonal.index[-1]].infil/12)
                    else:
                        df.loc[k]['Iw_'+tracer] = bar_d + (df.loc[k][tracer] - bar_d)* df.loc[k]['infil']/(df_seasonal.loc[df_seasonal.index[-1]].infil)
                    
                else:
                    print('cannot find annual infiltration sum, abort!', k)
                    df.loc[k]['Iw_'+tracer] = 0.
    return df

def alpha_deuterium(T):
    """a bit out of place here and should eventially be moved... alpha_d = alpha_dueterium(T)
        returns the isotopic fraction coefficient alpha defined from alpha = r_w/r_g
        where r_w is the isotopic ratio in the water phase and r_g is isotopic ratio
        in the gas phase.  The fraction factor is a function of T the temperature in Celcius
        and is interpolated from data from:
        J. Horita and D.J. Wesolowski, 1994, Liquid-vapour fractionation of oxygen and hydrogen isotopes of water
        from the freezing to the critical temperature, Geochimica et Cosmochimica Acta, 58, 3425-3437. """
    #table of dueterium partition coeficients from truesdell 1977
    #ln_alpha = array([106.,81.5,61.3,46.4,36.1,27.8,21.5,16.3,11.7,7.4,3.5,0.1,-2.2,-3.6,-4.,-3.4,-2.2,-1.3,-0.5]);
    #T_alpha = arange(0,361,20);
    #ln_alpha_T = interpolate.interpolate.spline(T_alpha,ln_alpha,T);
    T = T+273.15;
    # from Horita 1994
    alpha_horita_i =  ((1158.8*T**3)/10**9) - ((1640.1*T**2)/10**6) + ((794.84*T)/10**3) -161.04 + (2.9992e9/T**3);
    alpha = exp(alpha_horita_i/10**3);
    return alpha;

def alpha_oxygen(T):
    """a bit out of place here and should eventially be moved... alpha_d = alpha_dueterium(T)
        returns the isotopic fraction coefficient alpha defined from alpha = r_w/r_g
        where r_w is the isotopic ratio in the water phase and r_g is isotopic ratio
        in the gas phase.  The fraction factor is a function of T the temperature in Celcius
        and is interpolated from data from:
        J. Horita and D.J. Wesolowski, 1994, Liquid-vapour fractionation of oxygen and hydrogen isotopes of water
        from the freezing to the critical temperature, Geochimica et Cosmochimica Acta, 58, 3425-3437."""
    T = T+273.15;
    ln_alpha_i = 1137./(T**2) - 0.4156/T - 0.00207;
    alpha = exp(ln_alpha_i);
    return alpha;

def calc_alpha(df,gw_mean_conc,tracer='H2'):
    '''calculate alpha = alpha_s/alpla_w where alpha_i=I/P is the runoff ratio, using equation 2.15 from the IAEA book.  Needs the mean gw concentration measured in an area near the GNIP site. Tracer should be either 'H2' or 'O18'. For now assumes Northern hemisphere, and 2 6 month seasons. Summer march-september, winter october-feburary '''

    #add the precip weighted isotope P_i*d_i
    df.loc[:,'pw_tracer']=df[tracer]*df['Precipitation']
    #only use full "water years"
    prng = pd.period_range(datetime.datetime(df.index[0].year-1,3,1),datetime.datetime(df.index[-1].year+1,3,1),freq='12 M')
    df = df.groupby(prng.asof).filter(lambda x: len(x) == 12)
    #sum d*P and P for summer and winter
    SPd_s = 0.
    SPd_w = 0.
    SP_s = 0.
    SP_w = 0.
    for k in df.index:
        if k.month >= 3 and k.month <9:
            SPd_s = SPd_s + df.loc[k].pw_tracer
            SP_s = SP_s + df.loc[k].Precipitation
        else:  #winter
            SPd_w = SPd_w + df.loc[k].pw_tracer
            SP_w = SP_w + df.loc[k].Precipitation

    alpha = (SPd_w - gw_mean_conc*SP_w)/(gw_mean_conc*SP_s - SPd_s)
    return alpha



