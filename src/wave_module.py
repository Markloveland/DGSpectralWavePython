import numpy as np

#allocating variables
def init():
    global pi
    pi=np.pi
    global g
    g=9.81
    global deg_2_rad
    deg_2_rad=np.pi/180
    global rad_2_deg
    rad_2_deg=180/np.pi
    global dep_min
    dep_min=5.0000001*10**(-2)
    global run_ID
    run_ID="Empty"
    global cx_min, cy_min, csig_min, cth_min, cx_max, cy_max, csig_max, cth_max
    cx_min=float(99999)
    cy_min=float(99999)
    csig_min=float(99999)
    cth_min=float(99999)
    cx_max=float(-99999)
    cy_max=float(-99999)
    csig_max=float(-99999)
    cth_max=float(-99999)
    
    #Flags
    global stat_flag,div_sigma_flag,div_theta_flag,div_case_flag,theta_odd_flag,slope_flag
  
    stat_flag=99 #stationary or not
    div_sigma_flag=99 #if 0 the d csigmaN/dsigma =0
    div_theta_flag=99 #same as above but for theta
    div_case_flag=99 
    '''
    !DivCaseFlag = 1:DivSigmaFlag = 1, DivThetaFlag = 1
    !DivCaseFlag = 2:DivSigmaFlag = 1, DivThetaFlag = 0 
    !DivCaseFlag = 3:DivSigmaFlag = 0, DivThetaFlag = 1
    !DivCaseFlag = 0:DivSigmaFlag = 0, DivThetaFlag = 0
    '''
    theta_odd_flag=99 #0 is no change, 1 changes it to not coincide with grid
    slope_flag=99 #slope limiter
    global NF_override
    NF_override=99
    global full_circ_flag
    full_circ_flag=99
    global flag_wet_dry_elem
    flag_wet_dry_elem=99
    global wet_dry_flag
    wet_dry_flag=np.empty(99)
    global source_flag
    source_flag=99
    global p_tail
    p_tail=99 #calculate power of tail
    global out_x,out_y,out_depth,out_HS,out_DIR,out_TM01,out_TM02 #output flags
    out_x =99
    out_y=99
    out_depth=99
    out_HS=99
    out_DIR=99
    out_TM01=99
    out_TM02=99
    global out_RTM01,out_RTP,out_HS_max,out_HS_min,out_HS_dif
    out_RTM01=99
    out_RTP=99
    out_HS_max=99
    out_HS_min=99
    out_HS_dif=99
    global  sout_x,sout_y,sout_depth,sout_HS,sout_DIR,sout_TM01,sout_TM02 #station output flags
    sout_x=99
    sout_y=99
    sout_depth=99
    sout_HS=99
    sout_DIR=99
    sout_TM01=99
    sout_TM02=99
    global sout_RTM01,sout_RTP,sout_HS_max,sout_HS_min,sout_HS_dif,sout_HS_old,sout_DIR_old,sout_depth_old
    sout_RTM01=99
    sout_RTP=99
    sout_HS_max=99
    sout_HS_min=99
    sout_HS_dif=99
    sout_HS_old=99
    sout_DIR_old=99
    sout_depth_old=99

    #Variables for geographic mesh
    global num_elem, num_nodes, num_edges, num_inter_edges
    num_elem=99
    num_nodes=99
    num_edges=99
    num_inter_edges=99
    # array of internal edge numbers
    global inter_edge_array
    inter_edge_array=np.empty(99)
    global num_absorb_edges, absorb_edge_array #boundary absorbing edges
    num_absorb_edges=99
    absorb_edge_array=np.empty(99)
    #Incoming wave edges per boundary edge
    global num_in_wave_edges, in_wave_edge_array
    num_invwave_edges=np.empty(99)
    in_wave_edge_array=np.zeros((9,9))
    global num_refl_edges #number of reflective edges
    num_refl_edges=99
    global mnndel #max number of elements per node
    mnndel=99
    global node_2_elem #node to element table
    node_2_elem=np.zeros((9,9))
    

    ##structures for geographic mesh
    #WILL BE DONE LATER






    #error tolerance for convergence
    global epsilon_r, epsilon_a, epsilon_c,perct_elem
    epsilon_r=99
    epsilon_a=99
    epsilon_c=99
    perct_elem=99

    #variables for spectral mesh
    global f_min, f_max, d_theta #min and max frequency in Hz, spacing in radians
    f_min=float(99)
    f_max=float(99)
    d_theta=float(99)
    global num_sig_elem,theta_min,theta_max
    num_sig_elem=99
    theta_min=float(99)
    theta_max=float(99)
    global theta,sigma
    theta=np.empty(9)
    sigma=np.empty(9)
    global num_theta_elem
    num_theta_elem=99

    #variables for boundary conditions
    #number of incoming boundary condition segments
    global num_in_bc_seg
    num_in_bc_seg=99
    #segment number for which incoming wave is specified
    global in_bc_seg
    in_bc_seg=np.empty(99)
    global type_bc #1 pierson-moskowitz spectrum, 2 = jonswap spectrum, 4 = gauss
    type_bc=np.empty(99)
    global HS_bc #significant wave height in bc
    HS_bc=np.empty(99)
    global pk_per_bc #peak period in seconds
    pk_per_bc=np.empty(99)
    global sig_fr_bc #width of gaussian freq spectrum as std in Hz
    sig_fr_bc=np.empty(99)
    global m_dir_bc #wave direction in degrees for input radians
    m_dir_bc=np.empty(99)
    global ms_bc #power in directional distribution
    ms_bc=np.empty(99)
    global gamma_bc
    gamma_bc=np.empty(99)
    global n_bc #(qp,theta.sigma,qporder,bcsegment Action density N at boundary
    n_bc=np.empty((9,9,9,9,9))
    global n_bc_alpha #(k, theta, sigma, bcsegment)
    n_bc_alpha=np.empty((9,9,9,9))

    #Variables for DG things
    global p,p_low,p_high #order in geographic spacw
    p=99
    p_low=99
    p_high=99
    global q,q_low,q_high # spectral space
    q=99
    q_low=99
    q_high=99
    #dof in geo and spectral space
    global dof_XY, dof_SP
    dof_XY=np.empty(99)
    dof_SP=np.empty(99)
    global inv_M_mat #inverse of M matrix
    inv_M_mat=np.empty((99,99))
    global p_elem #p of a given element geographic
    p_elem=np.empty(99)
    global q_elem #q of specific element (theta,sigma, geographic)
    q_elem=np.empty((9,9,9))
    global RHS #RHS values
    RHS=np.empty((9,9,9,9,9,9))
    global alpha_N #coefficients of N=alpha*phi
    alpha_N=np.empty((9,9,9,9,9,9))
    global p_count, q_adapt
    p_count=np.empty(9)
    q_adapt=np.empty(9)

    #variables for time step
    global dt #time step in seconds
    dt=float(99)
    global T_time #number of time steps, if stationary max # iterations
    T_time=99
    global rad_stress_time_step
    rad_stress_time_step=99

    #variables for RK
    global RK_stage, RK_order
    RK_stage=99
    RK_order=99
    global a_RK_SSP, b_RK_SSP
    a_RK_SSP=np.empty((9,9))
    b_RK_SSP=np.empty((9,9))
    global delta_RK
    delta_RK=np.empty(9)
    global time_RK #time at RK step
    time_RK=float(9)

    #variables for depth and current
    global dof_u, p_u, read_current_64, time_varying_current
    dof_u=99
    p_u=99
    read_current_64=99
    time_varying_current=99
    global time_d1, time_d2
    time_d1=float(99)
    time_d2=float(99)
    global depth_time_step, depth_63_ic, num_data_sets_eta
    depth_time_step=99
    depth_63_ic=99
    num_data_sets_eta=99
    global time_c1, time_c2
    time_c1=float(99)
    time_c2=float(99)
    global current_time_step, current_64_ic, nof_data_sets_u
    current_time_step=99
    current_64_ic
    nof_data_sets_u=99
    global dof_d,p_d,read_eta_63,time_varying_depth
    dof_d=99
    p_d=99
    read_eta_63=99
    time_varying_depth=99
    global d_coef #depth coeff.
    d_coef=np.empty((9,9))
    global bath_coef
    bath_coef=np.empty((9,9))
    global ux_coef,uy_coef
    ux_coef=np.empty((9,9))
    uy_coef=np.wmpty((9,9))
    global aux,bux,cux,auy,buy,cuy #current in physical coord ux(x)-aux*x+bux*y+cux
    aux=float(99)
    bux=float(99)
    cux=float(99)
    auy=float(99)
    buy=float(99)
    cuy=float(99)
    global qx_coef,qy_coef
    qx_coef=np.empty((9,9))
    qy_coef=np.empty((9,9))

    #variables for quadrature points
    ##WILL BE COMPLETED LATER



    #variables for stations foroutput
    global num_sta,x_sta,y_sta,elem_sta,phi_sta
    num_sta = 99
    x_sta = np.empty(99)
    y_sta = np.empty(99)
    elem_sta=np.empty(99)
    phi_sta=np.empty((9,9))
    


def wave_alloc_boundary():
    global in_bc_seg
    global num_in_bc_seg
    in_bc_seg=np.empty(num_in_bc_seg)


