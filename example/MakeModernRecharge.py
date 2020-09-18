"""
Created on Mon Apr 30 15:12:45 2012
MakeConstraintList.py  writes a constraint list for PFloTran time series given a pandas time series dataset.
@author: wpgardn
"""
import pdb
from get_atm_conc import *

ddff = concat_resampled_time_series(C_t_gas,C_t_isotopes)
ddff_aq = convert2aqueous(ddff,25.0,addHe4=False,addKr81=False);
ddff_fin = convert2molality(ddff_aq);
#pdb.set_trace();
#ddff_fin.plot(logy=True);
#show();
make_modern_recharge(ddff_fin,model_time_0=ddff_fin.index[0])
