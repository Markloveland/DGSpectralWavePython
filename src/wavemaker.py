import wave_module
#import wave_mod_read
#import wave_messenger
#import wave_source
#import wave_source_quadruplets
#from wave_source_wind import wave_wind_initial
from wave_read_input import wave_ReadInput
#Already allocated in wave_module, just shown here for comparison sake with Jessica code
'''
wave_module.cx_min=float(99999)
wave_module.cy_min=float(99999)
wave_module.csig_min=float(99999)
wave_module.cth_min=float(99999)
wave_module.cx_max=float(-99999)
wave_module.cy_max=float(-99999)
wave_module.csig_max=float(-99999)
wave_module.cth_max=float(-99999)
'''

#read input parameter file
wave_ReadInput()
