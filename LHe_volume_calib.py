# Python code to plot the data from HELIX magnet test. Reads in .csv file downloaded form the HELIX wiki and plots versus time.  
# New code is to spline the levels on the near sensor of the magnet during the test to then take the derivative and plot as a function of time. 
#This is compared to the flowmeter data which is offset in time and converted to mass flow.  
# Author Keith McBride 11/20/18

import matplotlib.pyplot as plt
import numpy
from scipy import interpolate


##################################################################section on just splines of near level sensor
def load_levels_near(seq):
    datanear=numpy.genfromtxt('lvlSensorNear.csv', dtype=float, delimiter=',', names=True)
    timenear = datanear['time']
    lvlnear=datanear['StackSideLevelcm']
# convert unix time to seconds after tests started
    a=numpy.empty(len(timenear))
    a.fill(timenear[0])
    timenear=timenear-a
# convert to minutes after test started
    timenear= numpy.true_divide(timenear,60.0)
    data_list=numpy.array([timenear,lvlnear])
    return data_list

def spline_levels_near(seq):
# get the data
    data_list=load_levels_near(1)
# spline section
    #interpolate the data using spline of order 3 (cubic by default)
    tck, fp, ier, msg= interpolate.splrep(data_list[0,:], data_list[1,:], s=0, full_output=True)
    #new array for x-axis values to evaluate the calculated spline at
    time_spline=numpy.arange(0, data_list[0,-1], 0.1) #last argument can be changed for sampling time to better resolution
    #get the yvalues and the derivatives from the found spline at the xvalues from the new array above. 
    lvl_spline=interpolate.splev(time_spline, tck, der=0)
    lvl_der_spline=interpolate.splev(time_spline, tck, der=1)
#plot the spline
    fig=plt.figure(figsize=(10, 8), dpi=800)
    plt.scatter(time_spline,lvl_spline,c='r',marker='s', label='Near splined values')
    plt.scatter(data_list[0,:], data_list[1,:], c='b', marker='s', label='Near Sensor data')
    plt.legend(loc='upper right');
    plt.title('Magnet Thermal Test, near and splined near sensor levels')
    plt.ylabel('Level (cm)')
    plt.xlabel('Time (minutes)')
    fig.savefig('magnet_lvl_data_and_spline_vs_time.png')
    print (ier)
    print ('fp' ,fp)
    der_splined=numpy.array([time_spline,lvl_der_spline])
    return der_splined
#
def spline_der_levels_near(seq):
    der=spline_levels_near(1)
    fig=plt.figure(figsize=(10, 8), dpi=800)
    plt.scatter(der[0,:], der[1,:], c='b', marker='s', label='Near Sensor derivative')
    plt.title('Magnet Thermal Test, splined near sensor level derivative')
    plt.ylabel('Level derivative (cm per minute)')
    plt.xlabel('Time (mins)')
    fig.savefig('magnet_lvl_spline_der_vs_time.png')


################################################################### section below separate from above
def load_flows_near_lvl__and_der(seq):
# load lvlnear data
    datanear=numpy.genfromtxt('lvlSensorNear.csv', dtype=float, delimiter=',', names=True)
    timenear = datanear['time']
    lvlnear=datanear['StackSideLevelcm']
#load flow data
    data=numpy.genfromtxt('WhisperData_Cleaned.csv', dtype=float, delimiter=',', names=True)
    timeflow = data['time_s']
    pressure = data['pressure_PSI']
    temp = data['temperatureC']
    flow = data['volume_flow_LPM']
#    mass = data['mass_flow']
    #calculate the mass flow from the volume flow
      #convert temp to kelvin
    kelvin=numpy.empty(len(temp))
    kelvin.fill(273.15)
    temp=temp+kelvin
    masscalc = 4.0/(1000.0*1.20675)*numpy.true_divide(pressure,temp) #molar mass * 1000 * Pressure/(R[PSIliters/(kelvinmol)]*Temp)
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
#    data_list=numpy.array([timeflow,pressure,temp,flow,mass,masscalc])
    data_flow=numpy.array([timeflow,masscalc])
    data_near=numpy.array([timenear,lvlnear])
# just spline it here
    #interpolate the data using spline of order 3 (cubic by default)
    tck, fp, ier, msg= interpolate.splrep(data_near[0,:], data_near[1,:], s=0, full_output=True)
    #new array for x-axis values to evaluate the calculated spline at
    time_spline=numpy.arange(0, data_near[0,-1], 0.1) #last argument can be changed for sampling time to better resolution
    #get the yvalues and the derivatives from the found spline at the xvalues from the new array above. 
#    lvl_spline=interpolate.splev(time_spline, tck, der=0)
    lvl_der_spline=interpolate.splev(time_spline, tck, der=1)
    data_der=numpy.array([time_spline,lvl_der_spline])
   # end spline
    return data_flow, data_near, data_der


def plot_mass_flow_level_der(seq):
    print('starting')
#    der=spline_levels_near(1) # this loads the 2-d array of [0,:]= time (zeroed to near sensor start and in minutes) [1,:]= derivative of spline of level near sensor (cmpm)
    print('load everything')
#    flows=load_flows_zero_to_near_lvl_sensor(1)  # 5-d array with [0,:]=time (not zeroed to flow start and in minutes) [1,:]= pressure (PSI) [2,:]= temperature (Kelvin) [3,:]=volume flow (lpm) [4,:]= mass flow $
    flows, datanear, der=load_flows_near_lvl__and_der(1)
#    print('loaded flows')
#    datanear=load_levels_near(1) #[0,:]= time near from data taken (mins) [1,:]=lvl near cm reading (cm)
    print('about to make figure')
    # set flows time to zeroed 
    fig=plt.figure(figsize=(10, 8), dpi=800,)
    plt.title('Magnet Thermal Test, splined near sensor level derivative')
    ax_flow=plt.subplot(311)
    ax_flow.scatter(flows[0,:],flows[1,:],s=4)
    ax_flow.set_ylabel('Calculated massflow (kg per minute)')
    ax_flow.set_xlim([0,18000])
    ax_lvl=plt.subplot(312)
    ax_lvl.scatter(datanear[0,:],datanear[1,:],s=4)
    ax_lvl.set_xlim([0,18000])
    ax_lvl.set_ylabel('Level Near (cm)')
    ax_near=plt.subplot(313)
    ax_near.scatter(der[0,:],der[1,:],s=4)
    ax_near.set_xlim([0,18000])
    ax_near.set_ylabel('Level Near derivative(cm per min)')
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


def plot_flow_lvl_der(seq):
    plot_mass_flow_level_der(1)


############################################################section on using spline and mass flow for one plot


def find_section_of_lvl_sensor(seq):
    data_near= load_levels_near(1)
    #locate the element of data_near in which the level jumps?
    # in minutes it is after 2000 and before 4000. the jump occurs when the value is above 100.
    i=0
    times=[]
    lvls=[]
    while i< len(data_near[0,:]):
       if data_near[0,i]>2000 and data_near[0,i] <4000 and data_near[1,i]<100:
          times.append(data_near[0,i])
          lvls.append(data_near[1,i])
       i+=1
    data_splice=numpy.array([times,lvls])
    return data_splice

def plot_found_section(seq):
    section= find_section_of_lvl_sensor(1)
    fig=plt.figure(figsize=(10, 8), dpi=800)
    plt.scatter(section[0,:], section[1,:], c='b', marker='s', label='Near Sensor')
    plt.title('Magnet Thermal Test, near sensor section only')
    plt.ylabel('Level (cm)')
    plt.xlabel('Time (mins)')
    fig.savefig('magnet_lvl_sensor_vs_time_section.png')

def spline_section(seq):
    section= find_section_of_lvl_sensor(1)
# just spline it here
    #interpolate the data using spline of order 3 (cubic by default)
    tck, fp, ier, msg= interpolate.splrep(section[0,:], section[1,:], s=0, full_output=True)
    #new array for x-axis values to evaluate the calculated spline at
    time_spline=numpy.arange(section[0,0], section[0,-1], 0.2) #last argument can be changed for sampling time to better resolution
    #get the yvalues and the derivatives from the found spline at the xvalues from the new array above. 
# lvl_spline=interpolate.splev(time_spline, tck, der=0)
    lvl_der_spline=interpolate.splev(time_spline, tck, der=1)
    data_der=numpy.array([time_spline,lvl_der_spline])
    return data_der

def plot_spline_section(seq):
    der_spline=spline_section(1)
    fig=plt.figure(figsize=(10, 8), dpi=800)
    plt.scatter(der_spline[0,:], der_spline[1,:], c='b', marker='s', label='Near Sensor spline derivative')
    plt.title('Magnet Thermal Test, near sensor section only')
    plt.ylabel('Level derivative (cm per minute)')
    plt.xlabel('Time (mins)')
    fig.savefig('magnet_lvl_sensor_der_vs_time_section.png')

def find_section_flowmeter(seq):
    data_near= load_levels_near(1)
    #locate the element of data_near in which the level jumps?
    # in minutes it is after 2000 and before 4000. the jump occurs when the value is above 100.
    i=0
    times=[]
    lvls=[]
    while i< len(data_near[0,:]):
       if data_near[0,i]>2000 and data_near[0,i] <4000 and data_near[1,i]<100:
          times.append(data_near[0,i])
          lvls.append(data_near[1,i])
       i+=1
    data_splice=numpy.array([times,lvls])
    return data_splice


def plot_flowmeter_section(seq):
    
