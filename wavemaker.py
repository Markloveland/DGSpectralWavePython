import src.wave_module as wm
#import wave_mod_read
#import wave_messenger
#import wave_source
#import wave_source_quadruplets
#from wave_source_wind import wave_wind_ini
import src.wave_read_input as wr
import src.wave_test as tst
#Already allocated in wave_module, just shown here for comparison sake with Jessica code
'''
wm.cx_min=float(99999)
wm.cy_min=float(99999)
wm.csig_min=float(99999)
wm.cth_min=float(99999)
wm.cx_max=float(-99999)
wm.cy_max=float(-99999)
wm.csig_max=float(-99999)
wm.cth_max=float(-99999)
'''

#read input parameter file
wm.init()
wr.wave_ReadInput()

