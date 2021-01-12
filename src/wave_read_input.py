import src.wave_module as wm
#import wave_mod_read
#imporrt wave_messenger


def wave_ReadInput():
    #read input file called input_filenames where each line has the input file names
    input_filenames=open('input_filenames','r')
    #first line is name of parameter input file
    inputfile=input_filenames.readline().strip()
    #second line is grid in format of ADCIRC fort.14
    inputfile_geog=input_filenames.readline().strip()
    #third line is flag for if there are stations to record at or not
    is_sta=int(input_filenames.readline())
    if is_sta == 1:
        #read name of station input file name
        inputfile_sta=input_filenames.readline()
        wm.sout_x, wm.sout_y, wm.sout_depth, wm.sout_HS, wm.sout_DIR, wm.sout_TM01, wm.sout_TM02, wm.sout_RTM01, \
                wm.sout_RTP, wm.sout_HS_max, wm.sout_HS_min, wm.sout_HS_dif, wm.sout_HS_old, wm.sout_DIR_old, \
                wm.sout_depth_old = list(map(int, input_filenames.readline().split(',')))
    else:
        wm.num_sta, wm.sout_x, wm.sout_y, wm.sout_depth, wm.sout_HS, wm.sout_DIR, wm.sout_TM01, wm.sout_TM02, wm.sout_RTM01, \
                wm.sout_RTP, wm.sout_HS_max, wm.sout_HS_min, wm.sout_HS_dif, wm.sout_HS_old, wm.sout_DIR_old, \
                wm.sout_depth_old = [0]*16
    input_filenames.close()
    

    #now read parameter input file
    fin=open(inputfile,'r')
    #string with information about run
    wm.run_ID=fin.readline()
    #0 or 1 for stationary run
    wm.stat_flag=int(fin.readline().split('!')[0])
    if wm.stat_flag == 1:
        #percent of elements needed for convergence
        wm.perct_elem = float(fin.readline().split('!')[0])
        #error tolerance for relative, absolute, curvative error (HS)
        wm.epsilon_r, wm.epsilon_a, wm.epsilon_c = list(map(float, fin.readline().split('!')[0].split(',')))
    #read polynomial order in x p, p-low, p-high
    wm.p, wm.p_low, wm.p_high = list(map(int,fin.readline().split('!')[0].split(',')))
    #read polynomial order in spectral space
    wm.q, wm.q_low, wm.q_high = list(map(int,fin.readline().split('!')[0].split(',')))
    #read RK order, stage
    wm.RK_order, wm.RK_stage = list(map(int,fin.readline().split('!')[0].split(',')))
    #read slope limiter flag
    wm.slope_flag = int(fin.readline().split('!')[0])
    #time step in seconds
    wm.dt = float(fin.readline().split('!')[0])
    #Total number of time steps
    wm.T_time = int(fin.readline().split('!')[0])
    #Flag for calculating sigma derivative
    wm.div_sigma_flag = int(fin.readline().split('!')[0])
    #Flag for calculating theta derivative
    wm.div_theta_flag = int(fin.readline().split('!')[0])
    #Flag for changing direction in theta ???
    wm.theta_odd_flag = int(fin.readline().split('!')[0])
    #Non fatal over ride
    wm.NF_override = int(fin.readline().split('!')[0])
    #Wet Dry Flag
    wm.flag_wet_dry_elem = int(fin.readline().split('!')[0])
    #degrees of freedom for depth
    wm.p_d = int(fin.readline().split('!')[0])
    #read eta 63 file
    wm.read_eta_63 = int(fin.readline().split('!')[0])
    if wm.read_eta_63 == 1:
        wm.time_varying_depth, wm.depth_time_step, wm.depth_63_ic = list(map(int, fin.readline().split('!')[0].split(',')))
    #degree of freedom for current
    wm.p_u = int(fin.readline().split('!')[0])
    #1 = read 64 file for current, 0 = use aux, bux, etc (current in physical coord ux(x)=aux*x+bux*y+cux
    wm.read_current_64 = int(fin.readline().split('!')[0])
    if wm.read_current_64 == 0:
        wm.aux, wm.bux, wm.cux, wm.auy, wm.buy, wm.cuy =  list(map(float, fin.readline().split('!')[0].split(',')))
    elif wm.read_current_64 == 1:
        wm.time_varying_current, wm.current_time_step, wm.current_64_ic = list(map(int, fin.readline().split('!')[0].split(','))) 
    #info about spectral domain, calculate full circle or just sector
    wm.full_circ_flag = int(fin.readline().split('!')[0])
    if wm.full_circ_flag == 1:
        wm.num_theta_elem = int(fin.readline().split('!')[0])
        wm.theta_min=float(0)
        wm.theta_max=float(360)
    else:
        wm.num_theta_elem, wm.theta_min, wm.theta_max = list(map(float, fin.readline().split('!')[0].split(',')))
        wm.num_theta_elem = int(wm.num_theta_elem)
    wm.num_sig_elem,wm.f_min,wm.f_max= list(map(float, fin.readline().split('!')[0].split(',')))
    wm.num_sig_elem=int(wm.num_sig_elem)
    wm.out_x,wm.out_y,wm.out_depth,wm.out_HS,wm.out_DIR,wm.out_TM01,wm.out_TM02, \
            wm.out_RTM01,wm.out_RTP,wm.out_HS_max,wm.out_HS_min,wm.out_HS_dif= \
            list(map(int, fin.readline().split('!')[0].split(',')))
    wm.rad_stress_time_step=int(fin.readline().split('!')[0])
    wm.num_in_bc_seg=int(fin.readline().split('!')[0])
    wm.wave_alloc_boundary()
    print(wm.in_bc_seg)


    

    
