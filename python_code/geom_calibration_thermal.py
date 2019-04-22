# Python code to plot the data from HELIX magnet test. Reads in .csv file downloaded form the HELIX wiki and plots versus time the quantities.
# New code is to spline the levels on the magnet during the test to then take the derivative and plot as a function of time.
# Update including taking the ratio of the flow per time to the derivative of the near level sensor to get a cross-sectional area of the dewar as a function of the level sensor.
# Author Keith McBride 10/16/18

import matplotlib.pyplot as plt
import numpy
import paths
from scipy import interpolate
from matplotlib import rcParams
rcParams['mathtext.default'] = 'regular' #might be unneccessary
import objgraph


def load_flows_zero_to_near_lvl_sensor(seq):
# load lvlnear data
    datanear=numpy.genfromtxt(paths.data_path+'/lvlSensorNear.csv', dtype=float, delimiter=',', names=True)
    timenear = datanear['time']
    lvlnear=datanear['StackSideLevelcm']
#load flow data
    data=numpy.genfromtxt(paths.data_path+'/WhisperData_Cleaned.csv', dtype=float, delimiter=',', names=True)
    timeflow = data['time_s']
    pressure = data['pressure_PSI']
    temp = data['temperatureC']
    flow = data['volume_flow_LPM']
    mass = data['mass_flow']
    #calculate the mass flow from the volume flow
      #convert temp to kelvin
    kelvin=numpy.empty(len(temp))
    kelvin.fill(273.15)
    temp=temp+kelvin
    masscalc = 4.0/(1.20675)*numpy.true_divide(pressure,temp) #molar mass * Pressure/(R[PSIliters/(kelvinmol)]*Temp) has units of grams
    masscalc=numpy.multiply(masscalc,flow)
# convert unix time to seconds after tests started for all arrays of time, each array of time has different lengths since data was poorly recorded.
   # make array of same length as each time array
    a=numpy.empty(len(timeflow))
    b=numpy.empty(len(timenear))
   # fill each array of time with entry of earliest unix time (timenear[0])
    a.fill(timenear[0])
    b.fill(timenear[0])
  # subtract each start time array (a,b,c) from the corresponding time array.
    timeflow=timeflow-a
    timenear=timenear-b
# convert to minutes after test started
    timeflow= numpy.true_divide(timeflow,60.0)
    timenear= numpy.true_divide(timenear,60.0)
    data_list=numpy.array([timeflow,pressure,temp,flow,mass,masscalc])
    data_near=numpy.array([timenear,lvlnear])
    return data_list, data_near

def plot_mass_flow_and_level_der(seq):
    print('starting')
    der=spline_levels_near(1) # this loads the 2-d array of [0,:]= time (zeroed to near sensor start and in minutes) [1,:]= derivative of spline of level near sensor (cmpm)
    print('laoded derivative')
    flows=load_flows_zero_to_near_lvl_sensor(1)  # 5-d array with [0,:]=time (not zeroed to flow start and in minutes) [1,:]= pressure (PSI) [2,:]= temperature (Kelvin) [3,:]=volume flow (lpm) [4,:]= mass flow (? unknown units) [5,:]=mass flow calc from ideal gas law(kg per minute)
    print('loaded flows')
    datanear=load_levels_near(1) #[0,:]= time near from data taken (mins) [1,:]=lvl near cm reading (cm)
    print('abotu to make figure')
    # set flows time to zeroed
    fig=plt.figure(figsize=(10, 8), dpi=800,)
    plt.title('Magnet Thermal Test, splined near sensor level derivative')
    ax_flow=plt.subplot(311)
    ax_flow.scatter(flows[0,:],flows[5,:],s=4)
    ax_flow.set_ylabel('Calculated massflow (kg per minute)')
    ax_flow.set_xlim([0,18000])
    ax_near=plt.subplot(312)
    ax_near.scatter(der[0,:],der[1,:],s=4)
    ax_near.set_xlim([0,18000])
    ax_near.set_ylabel('Level Near derivative(cm per min)')
    ax_lvl=plt.subplot(313)
    ax_lvl.scatter(datanear[0,:],datanear[1,:],s=4)
    ax_lvl.set_xlim([0,18000])
    ax_lvl.set_ylabel('Level Near (cm)')
    ax_lvl.set_xlabel('Time (mins)')
    sum=0.0
    i=0
    while i<len(flows[0,:]):
      if i==len(flows[0,:])-1: diff=0
      else: diff=flows[0,i+1]-flows[0,i]
      sum+=diff
      i+=1
    print (sum/len(flows[0,:]))
#    plt.scatter(der[0,:], der[1,:], c='b', marker='s', label='Near Sensor derivative')
#    plt.scatter(der[0,:], der[1,:], c='r', marker='s', label='Near Sensor derivative')
    fig.savefig('magnet_mass_flow_calc_and_lvl_spline_der_vs_time.png')

######################for 11/28/18
# functions that will help find sections

# moving average for finding flows around level sensor values
def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def integrate_time_range(lvltime,flowtime,flow):
    flow_at_lvl_time=[]
#   integrate_timesA=[] #time units to integrate
    l=0
    find_times=[]
    while l<len(lvltime):
 #      print('lvl time before checking is ' , lvl_time_spliceA[l])
 #      print(l)
       #match up flows
       m=0
       if l==len(lvltime)-1:
          delta_t=lvltime[l]-lvltime[l-1]
       if l!=len(lvltime)-1:
          delta_t=lvltime[l+1]-lvltime[l]
       while m<len(flowtime):
          if m!=0 or m!=len(flowtime)-1: #not the end points of the array
             if flowtime[m]>lvltime[l]-(delta_t/2.0) and flowtime[m]<lvltime[l]+(delta_t/2.0):
#            if (time_spliceA[m]<lvl_time_spliceA[l] and time_spliceA[m+1]>lvl_time_spliceA[l]) or time_spliceA[m]==lvl_time_spliceA[l]:
                find_times.append(m)
#               print('time flow is' , time_spliceA[m])
#               print('lvl time is ' , lvl_time_spliceA[l])
#integrate from one time of lvlnear to the next to get an average flow.
#               average=(flows_spliceA[m-1]+flows_spliceA[m+1]+flows_spliceA[m])/3.0
          m+=1
       start=find_times[0]
       stop=find_times[-1]
       #print('start is ' , start)
       average=numpy.trapz(flow[start:stop], x=flowtime[start:stop])/(delta_t)
       #print (average)
       flow_at_lvl_time.append(average)
       find_times=[]
       l+=1
    return flow_at_lvl_time

def find_avg_errors(lvlA,lvlB,lvlC,ratA,ratB,ratC,num_bins):
    avg=[]
    errors=[]
    bins=[]
#stitch all the arrays together
    lvltots=[]
    rattots=[]
    l=0
    while l<len(lvlA):
       lvltots.append(lvlA[l])
       rattots.append(ratA[l])
       l+=1
    l=0
    while l<len(lvlB):
       lvltots.append(lvlB[l])
       rattots.append(ratB[l])
       l+=1
    l=0
    while l<len(lvlC):
       lvltots.append(lvlC[l])
       rattots.append(ratC[l])
       l+=1
#find the error and average in each bin
    #start with min and max of levels
    lvlmin=numpy.amin(lvltots)
    lvlmax=numpy.amax(lvltots)
    binwidth=(lvlmax-lvlmin)/num_bins
    print("{} is the binwidth in cm".format(binwidth))
    #outer loop for something about which bin
    lvliter=lvlmin-(binwidth/2.0)
    #array for storing the data points in the bin
    lvlbin_data=[]
    ratbin_data=[]
    while lvliter<lvlmax+(binwidth/2.0):
       #inner loop for checking all the elements
       l=0
       while l<len(lvltots):
          if lvltots[l]>=lvliter and lvltots[l]<=lvliter+(binwidth):
             lvlbin_data.append(lvltots[l]) # put that value in the array to be used for calculations
             ratbin_data.append(rattots[l])
             #find mean, reset values, declare values needed zero before this loop, and then calc standard deviation for those values somehow also
          l+=1
       print(len(ratbin_data))
#       if len(ratbin_data)==0:
#          ratbin_data.append(0)
#       avg.append(numpy.mean(ratbin_data))
#       errors.append(numpy.std(ratbin_data))
#       bins.append(lvliter)
       if len(ratbin_data)!=0:
          avg.append(numpy.mean(ratbin_data))
          errors.append(numpy.std(ratbin_data))
          bins.append(lvliter)
       lvliter+=binwidth
       lvlbin_data=[] # clear the bin data array
       ratbin_data=[]
    return bins, avg, errors

def plot_flows_multiple_subplots(seq):
####THIS FUNCTION AS IT IS PLOTS THE FLOW AND NEAR LEVEL SENSOR FOR ALL TIMES AND SHOWS THE SECTIONS CHOSEN FOR FURTHER ANALYSIS DUE TO EQ CONDITION
    flows, lvlnear=load_flows_zero_to_near_lvl_sensor(1)
###############plotter section############
###plot all###
    fig=plt.figure(figsize=(14, 14), dpi=800)
    if seq==1:
       fig.subplots_adjust(hspace=0)
       #need subplots, subplots(sharex=True)
       ax_flow=plt.subplot(211)
       plt.setp(ax_flow.get_xticklabels(), visible=False)
       ax_flow.scatter(flows[0,:],flows[3,:],color='k',s=30,marker='s')
       # add in the section lines
       xposition=[5138,6826.13,6900,8382.8,10000,16500]
       labels=[' A', ' B', 'section C']
       line_colors=['b','b','g','g','r','r']
       label_pos=[5300,7000,11800]
       label_style=['b','g','r']
       j=0
       i=0
       for xc in xposition:
          ax_flow.axvline(x=xc, color=line_colors[j], linestyle='--',linewidth=4) # this is for the lines to mark at what time notable things in the test occured.
          j+=1
       for xl in labels:
          ax_flow.text(label_pos[i], 40, labels[i], color=label_style[i],rotation=0, fontsize=32)
          i+=1
       ax_flow.set_ylabel('Flow (liters per minute)',fontsize=32)
       ax_flow.tick_params(axis='y',labelsize=28)
       ax_flow.set_xlim([0,18000])
       ax_flow.set_ylim([-9.9,60])
       ax_near=plt.subplot(212,sharex=ax_flow)
       plt.setp(ax_near.get_xticklabels(), rotation=20)
       ax_near.scatter(lvlnear[0,:],lvlnear[1,:],color='k',s=30, marker='s')
       j=0
       for xi in xposition:
          ax_near.axvline(x=xi, color=line_colors[j], linestyle='--',linewidth=4) # this is for the lines to mark at what time notable things in the test occured.
          j+=1
       ax_near.set_ylabel('Near Level (cm)',fontsize=32)
       ax_near.tick_params(axis='both', labelsize=28)
#       ax_near.tick_params(axis='x',tickdir=1)
       ax_near.set_xlim([0,17999])
       ax_near.set_ylim([0,139])
       ax_near.set_xlabel('Time (mins)',fontsize=32)
       fig.savefig(paths.images_path+'/Flow_and_lvl_vs_time_zeroed_together_sections_labelled.png')
       plt.clf()


def plot_sections_lvl_lvlsplined_and_derA(seq):
###load everything###
    flows, lvlnear=load_flows_zero_to_near_lvl_sensor(1)
###declare all necessary arrays###
    lvl_spliceA=[]
    lvl_time_spliceA=[]
    lvl_spliceB=[]
    lvl_time_spliceB=[]
    lvl_spliceC=[]
    lvl_time_spliceC=[]
    flows_spliceA=[]
    time_spliceA=[]
    flows_spliceB=[]
    time_spliceB=[]
    flows_spliceC=[]
    time_spliceC=[]
###splice flows into sections###
    k=0
    while k<len(flows[0,:]):
#        if ((flows[0,k]>5138 and flows[0,k]<5470) or (flows[0,k]>5490 and flows[0,k]<6826.13)) and flows[3,k]>20 and flows[3,k]<28: # amended this with getting rid of the data points in the window of 5470 and 5490
        if (flows[0,k]>5138 and flows[0,k]<5470 and flows[3,k]>20 and flows[3,k]<28) or (flows[0,k]>5490 and flows[0,k]<6826.13 and flows[3,k]>20 and flows[3,k]<28): # amended this with getting rid of the data points in the window of 5470 and 5490
           flows_spliceA.append(flows[5,k])
           time_spliceA.append(flows[0,k])
           #how to splice simultaneously the lvl sensor? its a smaller array so just go through the elements?
        if (flows[0,k]>6900 and flows[0,k]<7040) or (flows[0,k]<8382.8 and flows[0,k]>7060): # amended to avoid window of [7040,7060]
#        if flows[0,k]>6900 and flows[0,k]<8382.8:
           flows_spliceB.append(flows[5,k])
           time_spliceB.append(flows[0,k])
        if flows[0,k]>10000 and flows[0,k]<16500:
           flows_spliceC.append(flows[5,k])
           time_spliceC.append(flows[0,k])
        k+=1
###splice levels into sections###
    j=0
    while j<len(lvlnear[0,:]):
#        if lvlnear[0,j] >time_spliceA[0] and lvlnear[0,j]<time_spliceA[-1]:
        if (lvlnear[0,j] >time_spliceA[0] and lvlnear[0,j] <5470) or (lvlnear[0,j]<time_spliceA[-1] and lvlnear[0,j]>5490): #amended to avoid the window of 5470 and 5490
           lvl_spliceA.append(lvlnear[1,j])
           lvl_time_spliceA.append(lvlnear[0,j])
#        if lvlnear[0,j] >time_spliceB[0] and lvlnear[0,j]<time_spliceB[-1]:
        if (lvlnear[0,j]>time_spliceB[0] and lvlnear[0,j]<7040) or (lvlnear[0,j]>7060 and lvlnear[0,j]<time_spliceB[-1]): #amended to avoid the window of 7040 and 7060
           lvl_spliceB.append(lvlnear[1,j])
           lvl_time_spliceB.append(lvlnear[0,j])
        if lvlnear[0,j] >time_spliceC[0] and lvlnear[0,j]<time_spliceC[-1]:
           lvl_spliceC.append(lvlnear[1,j])
           lvl_time_spliceC.append(lvlnear[0,j])
        j+=1
#    print(time_spliceC[-1])
###spline sections###
    #sectionA spline
    tckA, fpA, ierA, msgA= interpolate.splrep(lvl_time_spliceA, lvl_spliceA, s=0, full_output=True)
    splineA_times_more=numpy.arange(lvl_time_spliceA[0],lvl_time_spliceA[-1],0.1)
    splineA=interpolate.splev(lvl_time_spliceA, tckA, der=0)
    splineA_more=interpolate.splev(splineA_times_more, tckA, der=0)
    splineA_der=interpolate.splev(lvl_time_spliceA, tckA, der=1)# derivative only at lvl times
    splineA_der_more=interpolate.splev(splineA_times_more, tckA, der=1)# derivative at spline times
    splineA_der_flowtimes=interpolate.splev(time_spliceA, tckA, der=1)# derivative at flow
    #sectionB spline
    tckB, fpB, ierB, msgB= interpolate.splrep(lvl_time_spliceB, lvl_spliceB, s=0, full_output=True)
    splineB_times_more=numpy.arange(lvl_time_spliceB[0],lvl_time_spliceB[-1],0.1)
    splineB=interpolate.splev(lvl_time_spliceB, tckB, der=0)
    splineB_more=interpolate.splev(splineB_times_more, tckB, der=0)
    splineB_der=interpolate.splev(lvl_time_spliceB, tckB, der=1)# derivative only at lvl times
    splineB_der_more=interpolate.splev(splineB_times_more, tckB, der=1)# derivative only at lvl times
    splineB_der_flowtimes=interpolate.splev(time_spliceB, tckB, der=1)# derivative at flow
    #sectionC spline
    tckC, fpC, ierC, msgC= interpolate.splrep(lvl_time_spliceC, lvl_spliceC, s=0, full_output=True)
    splineC_times_more=numpy.arange(lvl_time_spliceC[0],lvl_time_spliceC[-1],0.1)
    splineC=interpolate.splev(lvl_time_spliceC, tckC, der=0)
    splineC_more=interpolate.splev(splineC_times_more, tckC, der=0)
    splineC_der=interpolate.splev(lvl_time_spliceC, tckC, der=1) # derivative only at lvl times
    splineC_der_more=interpolate.splev(splineC_times_more, tckC, der=1)# derivative only at lvl times
    splineC_der_flowtimes=interpolate.splev(time_spliceC, tckC, der=1)# derivative at flow
###############plotter section############
###plot all###
#    fig=plt.figure(figsize=(16, 16), dpi=800)
    if seq==1:
#################################################for seq==1, plot all sections flows and lvls+splines.#############################################################
###plot sectionA###
    # plot flow
       fig, (ax_lvlA, ax_derA) = plt.subplots(figsize=(18,18), dpi=800, nrows=2,ncols=1,sharex='col')
       fig.subplots_adjust(hspace=0)
#       ((ax_lvlA,ax_derA),(ax_lvlB,ax_derB),(ax_lvlC,ax_derC))=plt.subplots(nrows=2,ncols=3,sharex='col')
#       plt.setp(ax_lvlA.get_xticklabels(), visible=False)
       ax_lvlA.scatter(splineA_times_more, splineA_more, c='k', s=20, label='Spline')
       ax_lvlA.scatter(lvl_time_spliceA, lvl_spliceA, c='b', s=80, label='Measurement')
       ax_lvlA.legend(loc='upper right',fontsize=28)
       ax_lvlA.tick_params(axis='y', labelsize=28)
       ax_lvlA.set_title('Section A',fontsize=32,color='b')
       ax_lvlA.set_ylim([76,95])
       ax_lvlA.set_ylabel('Level (cm)',fontsize=32)

       ax_derA.scatter(splineA_times_more, splineA_der_more, c='k', s=20, marker='s', label='Spline')
       ax_derA.scatter(lvl_time_spliceA, splineA_der, c='b', s=80, marker='s', label='Der. at data points')
       ax_derA.legend(loc='upper right', fontsize=28)
       ax_derA.ticklabel_format(style='sci',axis='y', scilimits=(-3,-3))
       ax_derA.yaxis.get_offset_text().set_visible(False)
       plt.text(4850, -0.002, "(1e-2)",fontsize=24)
#       ax_derA.yaxis.offsetText.set_fontsize(24)
#       ax_derA.yaxis.set_offset_position('right')
       ax_derA.set_ylim([-0.02,-0.001])
       ax_derA.set_xlim([5001,7000])
       ax_derA.tick_params(axis='both', labelsize=28)
#       ax_derA.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True))
       ax_derA.set_ylabel('Rate of level change (cm per min)',fontsize=32)
       ax_derA.set_xlabel('Time (mins)',fontsize=32)

       fig.savefig(paths.images_path+'/All_sections_lvl_lvlspline_and_lvlder_vs_time.png')
   #    fig.savefig('sectionA_lvl_vs_time.eps')
       plt.clf()
      #plot flow and level with subplots?
   #https://matplotlib.org/devdocs/gallery/subplots_axes_and_figures/shared_axis_demo.html#sphx-glr-gallery-subplots-axes-and-figures-shared-axis-demo-py
   #    ax_spline=plt.subplot(211,sharex


def load_all_arrays(seq):
###load everything###
    flows, lvlnear=load_flows_zero_to_near_lvl_sensor(1)
###declare all necessary arrays###
    lvl_spliceA=[]
    lvl_time_spliceA=[]
    lvl_spliceB=[]
    lvl_time_spliceB=[]
    lvl_spliceC=[]
    lvl_time_spliceC=[]
    flows_spliceA=[]
    time_spliceA=[]
    flows_spliceB=[]
    time_spliceB=[]
    flows_spliceC=[]
    time_spliceC=[]
###splice flows into sections###
    k=0
    while k<len(flows[0,:]):
#        if ((flows[0,k]>5138 and flows[0,k]<5470) or (flows[0,k]>5490 and flows[0,k]<6826.13)) and flows[3,k]>20 and flows[3,k]<28: # amended this with getting rid of the data points in the window of 5470 and 5490
        if (flows[0,k]>5138 and flows[0,k]<5470 and flows[3,k]>20 and flows[3,k]<28) or (flows[0,k]>5490 and flows[0,k]<6826.13 and flows[3,k]>20 and flows[3,k]<28): # amended this with getting rid of the data points in the window of 5470 and 5490
           flows_spliceA.append(flows[5,k])
           time_spliceA.append(flows[0,k])
           #how to splice simultaneously the lvl sensor? its a smaller array so just go through the elements?
        if (flows[0,k]>6900 and flows[0,k]<7040) or (flows[0,k]<8382.8 and flows[0,k]>7060): # amended to avoid window of [7040,7060]
#        if flows[0,k]>6900 and flows[0,k]<8382.8:
           flows_spliceB.append(flows[5,k])
           time_spliceB.append(flows[0,k])
        if flows[0,k]>10000 and flows[0,k]<16500:
           flows_spliceC.append(flows[5,k])
           time_spliceC.append(flows[0,k])
        k+=1
###splice levels into sections###
    j=0
    while j<len(lvlnear[0,:]):
#        if lvlnear[0,j] >time_spliceA[0] and lvlnear[0,j]<time_spliceA[-1]:
        if (lvlnear[0,j] >time_spliceA[0] and lvlnear[0,j] <5470) or (lvlnear[0,j]<time_spliceA[-1] and lvlnear[0,j]>5490): #amended to avoid the window of 5470 and 5490
           lvl_spliceA.append(lvlnear[1,j])
           lvl_time_spliceA.append(lvlnear[0,j])
#        if lvlnear[0,j] >time_spliceB[0] and lvlnear[0,j]<time_spliceB[-1]:
        if (lvlnear[0,j]>time_spliceB[0] and lvlnear[0,j]<7040) or (lvlnear[0,j]>7060 and lvlnear[0,j]<time_spliceB[-1]): #amended to avoid the window of 7040 and 7060
           lvl_spliceB.append(lvlnear[1,j])
           lvl_time_spliceB.append(lvlnear[0,j])
        if lvlnear[0,j] >time_spliceC[0] and lvlnear[0,j]<time_spliceC[-1]:
           lvl_spliceC.append(lvlnear[1,j])
           lvl_time_spliceC.append(lvlnear[0,j])
        j+=1
#    print(time_spliceC[-1])
###spline sections###
    #sectionA spline
    tckA, fpA, ierA, msgA= interpolate.splrep(lvl_time_spliceA, lvl_spliceA, s=0, full_output=True)
    splineA_times_more=numpy.arange(lvl_time_spliceA[0],lvl_time_spliceA[-1],0.1)
    splineA=interpolate.splev(lvl_time_spliceA, tckA, der=0)
    splineA_more=interpolate.splev(splineA_times_more, tckA, der=0)
    splineA_der=interpolate.splev(lvl_time_spliceA, tckA, der=1)# derivative only at lvl times
    splineA_der_more=interpolate.splev(splineA_times_more, tckA, der=1)# derivative at spline times
    splineA_der_flowtimes=interpolate.splev(time_spliceA, tckA, der=1)# derivative at flow
    #sectionB spline
    tckB, fpB, ierB, msgB= interpolate.splrep(lvl_time_spliceB, lvl_spliceB, s=0, full_output=True)
    splineB_times_more=numpy.arange(lvl_time_spliceB[0],lvl_time_spliceB[-1],0.1)
    splineB=interpolate.splev(lvl_time_spliceB, tckB, der=0)
    splineB_more=interpolate.splev(splineB_times_more, tckB, der=0)
    splineB_der=interpolate.splev(lvl_time_spliceB, tckB, der=1)# derivative only at lvl times
    splineB_der_more=interpolate.splev(splineB_times_more, tckB, der=1)# derivative only at lvl times
    splineB_der_flowtimes=interpolate.splev(time_spliceB, tckB, der=1)# derivative at flow
    #sectionC spline
    tckC, fpC, ierC, msgC= interpolate.splrep(lvl_time_spliceC, lvl_spliceC, s=0, full_output=True)
    splineC_times_more=numpy.arange(lvl_time_spliceC[0],lvl_time_spliceC[-1],0.1)
    splineC=interpolate.splev(lvl_time_spliceC, tckC, der=0)
    splineC_more=interpolate.splev(splineC_times_more, tckC, der=0)
    splineC_der=interpolate.splev(lvl_time_spliceC, tckC, der=1) # derivative only at lvl times
    splineC_der_more=interpolate.splev(splineC_times_more, tckC, der=1)# derivative only at lvl times
    splineC_der_flowtimes=interpolate.splev(time_spliceC, tckC, der=1)# derivative at flow

#####PUT ALL ARRAYS INTO ONE ARRAY######
    lvl_time_splice=[]
    lvl_splice=[]
    spline_times_more=[]
    spline_more=[]
    spline_der=[]
    spline_der_more=[]
    lvl_time_splice.append(lvl_time_spliceA)
    lvl_time_splice.append(lvl_time_spliceB)
    lvl_time_splice.append(lvl_time_spliceC)
    lvl_splice.append(lvl_spliceA)
    lvl_splice.append(lvl_spliceB)
    lvl_splice.append(lvl_spliceC)
    spline_times_more.append(splineA_times_more)
    spline_times_more.append(splineB_times_more)
    spline_times_more.append(splineC_times_more)
    spline_more.append(splineA_more)
    spline_more.append(splineB_more)
    spline_more.append(splineC_more)
    spline_der.append(splineA_der)
    spline_der.append(splineB_der)
    spline_der.append(splineC_der)
    spline_der_more.append(splineA_der_more)
    spline_der_more.append(splineB_der_more)
    spline_der_more.append(splineC_der_more)
#    all=[lvl_time_splice, lvl_splice, spline_times_more, spline_more, spline_der, spline_der_more]
    if (seq==1):
       return lvl_time_splice, lvl_splice, spline_times_more, spline_more, spline_der, spline_der_more
    else:
       time_splice=[]
       flows_splice=[]
       time_splice.append(time_spliceA)
       time_splice.append(time_spliceB)
       time_splice.append(time_spliceC)
       flows_splice.append(flows_spliceA)
       flows_splice.append(flows_spliceB)
       flows_splice.append(flows_spliceC)
       return lvl_time_splice, lvl_splice, spline_times_more, spline_more, spline_der, spline_der_more, flows_splice, time_splice
#    return all
def plot_sections_lvl_lvlsplined_and_der_all(seq):
    lvl_time_splice, lvl_splice, spline_times_more, spline_more, spline_der, spline_der_more= load_all_arrays(1)
#    all=load_all_arrays(1)
    objgraph.show_most_common_types()

###############plotter section############
###plot all###
#    fig=plt.figure(figsize=(16, 16), dpi=800)
    if seq==1:
#################################################for seq==1, plot all sections lvls+splines.#############################################################
###plot sections###
    # plot flow
       fig, (ax_lvl, ax_der) = plt.subplots(figsize=(24,10), dpi=1200, nrows=2,ncols=3,sharex='col')
       fig.subplots_adjust(hspace=0)
       objgraph.show_most_common_types()

#       ((ax_lvlA,ax_derA),(ax_lvlB,ax_derB),(ax_lvlC,ax_derC))=plt.subplots(nrows=2,ncols=3,sharex='col')
#       plt.setp(ax_lvlA.get_xticklabels(), visible=False)
       colors_list=['b','g','r']
       title_list=['Section A','Section B','Section C']
       y_lim_lvl=[[76,95],[65,77],[28,70]]
       y_lim_der=[[-0.025,0.006],[-0.025,0.006],[-0.025,0.006]]
       x_lim_der=[[4999,7000],[6802,8601],[9002,17001]]
       txt_offsetx=[-11268,-1675,7930]
       txt_offsety=[0.002,0.002,0.002]
       i=0
       while i < 3:
          print (len(spline_times_more[i]))
          ax_lvl[i].scatter(spline_times_more[i], spline_more[i], c='k', s=20, label='Spline')
          ax_lvl[i].scatter(lvl_time_splice[i], lvl_splice[i], c=colors_list[i], s=80, label='Measurement')
          ax_lvl[i].legend(loc='upper right',fontsize=14)
          ax_lvl[i].tick_params(axis='y', labelsize=18)
          ax_lvl[i].set_title(title_list[i],fontsize=24,color=colors_list[i])
          ax_lvl[i].set_ylim(y_lim_lvl[i])
          ax_lvl[i].set_ylabel('Level (cm)',fontsize=18)

          ax_der[i].scatter(spline_times_more[i], spline_der_more[i], c='k', s=20, marker='s', label='Spline')
          ax_der[i].scatter(lvl_time_splice[i], spline_der[i], c=colors_list[i], s=80, marker='s', label='Der. at data points')
          ax_der[i].legend(loc='lower right', fontsize=14)
          ax_der[i].ticklabel_format(style='sci',axis='y', scilimits=(-2,-2))
          plt.setp(ax_der[i].get_xticklabels(), rotation=20)
          ax_der[i].yaxis.get_offset_text().set_visible(False)
          plt.text(txt_offsetx[i],txt_offsety[i], "(1e-2)",fontsize=18)
#       ax_derA.yaxis.offsetText.set_fontsize(24)
#       ax_derA.yaxis.set_offset_position('right')
          ax_der[i].set_ylim(y_lim_der[i])
          ax_der[i].set_xlim(x_lim_der[i])
          ax_der[i].tick_params(axis='both', labelsize=18)
#       ax_derA.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True))
          ax_der[i].set_ylabel('Rate of change (cm per min)',fontsize=18)
          ax_der[i].set_xlabel('Time (mins)',fontsize=18)
          i+=1

       fig.savefig(paths.images_path+'/All_sections_lvl_lvlspline_and_lvlder_vs_time.png')
   #    fig.savefig('sectionA_lvl_vs_time.eps')
       plt.clf()
       plt.close('all')
      #plot flow and level with subplots?
   #https://matplotlib.org/devdocs/gallery/subplots_axes_and_figures/shared_axis_demo.html#sphx-glr-gallery-subplots-axes-and-figures-shared-axis-demo-py
   #    ax_spline=plt.subplot(211,sharex
def plot_ratio(seq):
    #Find the ratio of the flow and the level derivative to get the cross-sectional vs time measure of the dewar.
    lvl_time_splice, lvl_splice, spline_times_more, spline_more, spline_der, spline_der_more, flows_splice, time_splice= load_all_arrays(2)
       #calc the average at each time of lvl sensor data points for each section at once
    flow_at_lvl=[]
    ratio=[]
    j=0
    while j < 3:
       spline_der[j]=spline_der[j]*(-1.0)
       flow_integrated=integrate_time_range(lvl_time_splice[j], time_splice[j], flows_splice[j])
       flow_at_lvl.append(flow_integrated)
       print(len(flow_at_lvl[j]))
       print(len(lvl_time_splice[j]))
       ratio.append(numpy.true_divide(flow_integrated,spline_der[j]))
       ratio[j]=ratio[j]/1250.0 # density of the LHe is 0.125 g per cm cube. So divide by this and convert to meters squared
       # convert from grams per cm to m squared
       j+=1
    #now get rid of bad values (that derivative is not calculated correctly).
    j=0
    ratio_refined=[]
    lvl_refined=[]
    lvl_time_refined=[]
    while j < 3:
       n=0
       ratio_intermediate=[]
       lvl_intermediate=[]
       lvl_time_intermediate=[]
       while n<len(ratio[j]):
          if ratio[j][n]<1 and ratio[j][n]>0:
             ratio_intermediate.append(ratio[j][n])
             lvl_intermediate.append(lvl_splice[j][n])
             lvl_time_intermediate.append(lvl_time_splice[j][n])
          n+=1
       ratio_refined.append(ratio_intermediate)
       lvl_refined.append(lvl_intermediate)
       lvl_time_refined.append(lvl_time_intermediate)
       j+=1
    j=40
    bins=[]
    avg=[]
    errors=[]
    p_bins=[]
    plot_name=[]
    while j>=5:
       bins_i, avg_i, errors_i=find_avg_errors(lvl_refined[0],lvl_refined[1],lvl_refined[2],ratio_refined[0],ratio_refined[1],ratio_refined[2],j)
       bins.append(bins_i)
       avg.append(avg_i)
       errors.append(errors_i)
       p_bins.append(numpy.poly1d(numpy.polyfit(bins_i, avg_i, 0, full=False)))
       plot_name.append("{delta}{binwidth}".format(delta=r'$\Delta L=$',binwidth=bins_i[1]-bins_i[0]))
       j=j/2

    fntsz=24 #fontsize is specified
    width_of_line=3
    size_of_marker=6
    colors_list=['r','b','g','orange']
    y_lim=[0.01,0.99]
    x_axis_label=[53,-0.3]
    y_axis_label=[10,2.4]
    #RECONFIGURE INTO LOOP OVER AXES SUBPLOTS
    fig, ax_bins = plt.subplots(figsize=(14,14), dpi=1200, nrows=4,ncols=1,sharex='col')
    fig.subplots_adjust(hspace=0)
    j=0
    while j <= 3:
       ax_bins[j].errorbar(bins[j], avg[j], errors[j], linestyle='None', marker='^', elinewidth=size_of_marker, color=colors_list[j],label=plot_name[j])
       ax_bins[j].axhline(p_bins[j](50),color='black', lw=width_of_line,label='Constant Fit')
       ax_bins[j].legend(loc='upper left',fontsize=fntsz-4)
       ax_bins[j].set_ylim(y_lim)
       ax_bins[j].tick_params(axis='both', labelsize=fntsz-4)
       j=j+1
    plt.text(x_axis_label[0],x_axis_label[1], 'Level (cm)',fontsize=fntsz)
    plt.text(y_axis_label[0],y_axis_label[1], 'Cross-section (m$^2$)',fontsize=fntsz,rotation=90)
    fig.suptitle('Ratio vs Level, Various Binnings',fontsize=fntsz)
    filename="Ratio_vs_levels_binwidths_while_loop.png"
    fig.savefig(paths.images_path+'/'+filename)
##Nick's code
# generate plotting arrays for volume vs lvl plot. 
    plt.clf()
    cryo_data=numpy.genfromtxt(paths.data_path+'/vol_vs_lvl_cryo.csv', dtype=float, delimiter=',', names=True)
    lvl_cryo = cryo_data['lvl']
    vol_cryo = cryo_data['vol']
    x=numpy.arange(lvl_cryo[2],100,1)
    j=0
    fit=[]
    slope=[]
    while j <= 3:
       fit.append((p_bins[j](50)*10.0)*x+(vol_cryo[2]-(p_bins[j](50)*10.0*x[0])))
       slope.append(format(p_bins[j](50)*10.0,'4.3f'))
       j+=1
#    plt.xlim(30,100)
#    plt.ylim(0,300)
    nvol=[]
    nlvl=[]
    j = 0
    while j < len(lvl_cryo):
        if lvl_cryo[j]<=100 and lvl_cryo[j]>=30:
            nlvl.append(lvl_cryo[j])
            nvol.append(vol_cryo[j])
            print(nlvl[-1])
        j+=1
    x2 =numpy.arange(nlvl[0],nlvl[-1],1)
    y2=numpy.polyfit(nlvl,nvol,1)
    y2[1]=y2[1]+nvol[0]
    z2=numpy.polyval(y2,x2)
    slope.append(format(y2[0],'4.3f'))
#    plt.plot(x,fit[0],'-.',x,fit[1],':',x,fit[2],x,fit[3],'--',nlvl,nvol,'o',x2,z2,dashes=[6,2], lw=width_of_line)
    plt.plot(x,fit[0],'-.', lw=width_of_line)
    plt.plot(x,fit[1],':', lw=width_of_line)
    plt.plot(x,fit[2], lw=width_of_line)
    plt.plot(x,fit[3],'--', lw=width_of_line)
    plt.plot(nlvl,nvol,'o', markersize=size_of_marker)
    plt.plot(x2,z2,dashes=[6,2], lw=width_of_line)

#    plt.plot(x2,z2,c='red',dashes=[6,2])
    plt.xlim(29,100)
    plt.ylim(0,300)
    plt.legend(("m="+slope[0]+" "+plot_name[0],"m="+slope[1]+" "+plot_name[1],"m="+slope[2]+" "+plot_name[2],"m="+slope[3]+" "+plot_name[3],'Cryo Table',"m="+slope[4]+"  Fit Cryo table"),loc='upper left',prop={'size': fntsz-4})
    plt.title('Geometric calibration using gaseous flow',fontsize=fntsz)
    plt.ylabel('Volume (liters)',fontsize=fntsz)
    plt.xlabel('Sensor Level (cm)',fontsize=fntsz)
    plt.tick_params(labelsize=fntsz-4)
    filename="Volume_vs_levels_cryo_doc.png"
    fig.savefig(paths.images_path+'/'+filename)
