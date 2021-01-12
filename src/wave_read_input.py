#import wave_module
#import wave_mod_read
#imporrt wave_messenger




def wave_ReadInput():
    #read input file called input_filenames where each line has the input file names
    input_filenames=open('input_filenames','r')
    #first line is name of parameter input file
    inputfile=input_filenames.readline()
    #second line is grid in format of ADCIRC fort.14
    inputfile_geog=input_filenames.readline()
    #third line is flag for if there are stations to record at or not
    is_sta=int(input_filenames.readline())
    if is_sta == 1:
        #read name of station input file name
        inputfile_sta=input_filenames.readline()
        sout_x, sout_y, sout_depth, sout_HS, sout_DIR, sout_TM01, sout_TM02, sout_RTM01, \
                sout_RTP, sout_HSmax, sout_HSmin, sout_HSdif, sout_HSold, sout_DIRold, \
                sout_depthold = list(map(int, input_filenames.readline().split(',')))
    else:
        NumSta, sout_x, sout_y, sout_depth, sout_HS, sout_DIR, sout_TM01, sout_TM02, sout_RTM01, \
                sout_RTP, sout_HSmax, sout_HSmin, sout_HSdif, sout_HSold, sout_DIRold, \
                sout_depthold = [0]*16

    #now read parameter input file
    fin=open(inputfile,'r')
    #string with information about run
    run_ID=fin.readline()
    #0 or 1 for stationary run
    stat_flag=int(fin.readline())
    if stat_flag == 1:
        #percent of elements needed for convergence
        perct_elem = float(fin.readline())
        #error tolerance for relative, absolute, curvative error (HS)
        epsilon_r, epsilon_a, epsilon_c = list(map(float, fin.readline().split(',')))
    #read polynomial order in x p, p-low, p-high
    p, p_low, p_high = list(map(int,fin.readline().split(',')))
    #read polynomial order in spectral space
    q, q_low, q_high = list(map(int,fin.readline().split(',')))
    #read RK order, stage
    RK_order, RK_stage = list(map(int,fin.readline().split(',')))
    #read slope limiter flag
    slope_flag = int(fin.readline())
    #time step in seconds
    dt = float(fin.readline())
    

    
