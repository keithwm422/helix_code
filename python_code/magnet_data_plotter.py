# Python code to plot the data from HELIX magnet test. Reads in .csv file downloaded form the HELIX wiki and plots versus time the quantities.  
# New code is to spline the levels on the magnet during the test to then take the derivative and plot as a function of time. 
# Author Keith McBride 10/16/18

import matplotlib.pyplot as plt
import numpy
from scipy import interpolate
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

def load_file(seq):
    data=numpy.genfromtxt('WhisperData_Cleaned.csv', dtype=float, delimiter=',', names=True)
    time = data['time_s']
# convert unix time to seconds after tests started
    a=numpy.empty(len(time))
    a.fill(time[0])
    time=time-a    
# convert to hours after test started
    time= numpy.true_divide(time,3600.0)
    pressure = data['pressure_PSI']
    temp = data['temperatureC']
    flow = data['volume_flow_LPM']
    mass = data['mass_flow']
    #calculate the mass flow from the volume flow
      #convert temp to kelvin
    kelvin=numpy.empty(len(temp))
    kelvin.fill(273.15)
    temp=temp+kelvin
    masscalc = 4.0/(1000.0*1.20675)*numpy.true_divide(pressure,temp) #molar mass * 1000 * Pressure/(R[PSIliters/(kelvinmol)]*Temp)
    masscalc=numpy.multiply(masscalc,flow)
#Can add plot code here to load and save the flux plot for debugging later for the user. will be helpful for HELIX to have this universal with labels. 
#    data=numpy.array(data)
#    i=0
#    while i < len(data):
#      time.append(data[i][0])
#      pressure.append(data[i][1])
#      temp.append(data[i][2])
#      volume.append(data[i][3])
#      mass.append(data[i][4])
#      i+=1
    data_list=numpy.array([time,pressure,temp,flow,mass,masscalc])
    return data_list

def calc_average_density(seq):
    data_list=load_file(1)
    i=0
    tot=0
    sum=0
    while i < len(data_list[3,:]):
      if data_list[3,i]!=0:
        sum+=data_list[4,i]/data_list[3,i]
        tot+=1
      i+=1
    print('average density used is {}'.format(sum/tot)) 
#    fig=plt.figure(figsize=(10, 8), dpi=800,)
#    plt.scatter(data_list[0,:],density, s=4)
#    plt.title('Magnet whisper testing')
#    plt.ylabel('density(mass over liters)')
#    plt.show()
def integrate_mass(seq):
    data_list=load_file(1)
    #convert time to minutes instead of hours
    data_list[0,:]= 60.0*data_list[0,:]
    
    #integrate over time and report the value
    tot=numpy.trapz(data_list[5,:],x=data_list[0,:])
#    i=0
#    tot=0
#    while i<len(data_list[4,:]):
      
      
#      i+=1
    print('total mass flowed out is {}'.format(tot)) 

def integrate_volume_range(seq):
    data_list=load_file(1)
    fig=plt.figure(figsize=(10, 8), dpi=800,)
#    start = 2170 #found from grep -n "1510520760" WhisperCleaned_data.csv
#    end = 251968 # found from same as above but for time of 1510854954. also, because of 1st line in WhisperData_cleaned.csv being headers, subtract one from the line number. Then subtract one again for the accessing elements of an array (start at zero). 
    start = 251968 # found from same as above but for time of 1510520760. also, because of 1st line in WhisperData_cleaned.csv being headers, subtract one from the line number. Then subtract one again for the accessing elements of an array (start at zero). 
    end = len(data_list[0,:])-1 # end of recording flows 

    #convert time to minutes instead of hours
    data_list[0,:]= numpy.true_divide(data_list[0,:],0.016666666667)
    
    #integrate over time specified here and report the value
    tot=numpy.trapz(data_list[3,start:end],x=data_list[0,start:end])
    plt.scatter(data_list[0,start:end],data_list[3,start:end],s=4)
    plt.title('Magnet Mass flow Test time range')
    plt.ylabel('Volume (LPM)')
    plt.xlabel('Time (hours)')

    fig.savefig('volume_lpm_vs_time.png')

#    i=0
#    tot=0
#    while i<len(data_list[4,:]):
#      i+=1
#    print('total mass flowed out from Nov12th to Nov 16th is {}'.format(tot)) 
    print('total volume flowed out from Nov16th to Nov 21st is {}'.format(tot)) 
    tot_mass_const_density=tot*(0.15) #at this temp, density of helium gas in grams per liter
    print('total mass flowed out from Nov16th to Nov 21st is {} grams'.format(tot_mass_const_density))    


def integrate_mass_range(seq):
    data_list=load_file(1)
    fig=plt.figure(figsize=(10, 8), dpi=800,)
#    start = 2170 #found from grep -n "1510520760" WhisperCleaned_data.csv
#    end = 251968 # found from same as above but for time of 1510520760. also, because of 1st line in WhisperData_cleaned.csv being headers, subtract one from the line number. Then subtract one again for the accessing elements of an array (start at zero). 
    start = 251968 # found from same as above but for time of 1510520760. also, because of 1st line in WhisperData_cleaned.csv being headers, subtract one from the line number. Then subtract one again for the accessing elements of an array (start at zero). 
    end = len(data_list[0,:])-1 # found from same as above but for time of 1510520760. also, because of 1st line in WhisperData_cleaned.csv being headers, subtract one from the line number. Then subtract one again for the accessing elements of an array (start at zero). 

    #convert time to minutes instead of hours
    data_list[0,:]= numpy.true_divide(data_list[0,:],0.016666666667)
    
    #integrate over time specified here and report the value
    tot=numpy.trapz(data_list[4,start:end],x=data_list[0,start:end])
    plt.scatter(data_list[0,start:end],data_list[4,start:end],s=4)
    plt.title('Magnet Mass flow Test time range')
    plt.ylabel('Mass (MPM)')
    plt.xlabel('Time (hours)')

    fig.savefig('mass_mpm_vs_time.png')

#    i=0
#    tot=0
#    while i<len(data_list[4,:]):
#      i+=1
#    print('total mass flowed out from Nov12th to Nov 16th is {}'.format(tot)) 
    print('total mass flowed out from Nov16th to Nov 21st is {}'.format(tot)) 


def plot_flow_all(seq):
    data_list=load_file(1);
    time0=2170
    time1=251968
    time2=3289
#    plt.yscale('log')
    fig=plt.figure(figsize=(10, 8), dpi=800,)

    plt.subplot(511)
    plt.scatter(data_list[0,0:-1:60],data_list[1,0:-1:60], s=4)

# need to put these lines (but not the labels) on each subplot. Need to resize the text and add more lines for other notable things during data taking. 
#    xposition = [data_list[0,time0], data_list[0,time2], data_list[0,time1], data_list[0,-1]] # these will be values in hours.
#    labels=['Ramp start, FULL', 'start ramping ','Ramp end, FULL','end of data']
#    i=0
#    for xc in xposition:
#      plt.axvline(x=xc, color='k', linestyle='--') # this is for the lines to mark at what time notable things in the test occured.
#      plt.text(xc+1, 10, labels[i], rotation=90, fontsize=6)
#      i+=1

    plt.title('Magnet Thermal Test')
    plt.ylabel('Pressure (PSI)')
    plt.subplot(512)
    plt.scatter(data_list[0,0:-1:60],data_list[2,0:-1:60], s=4)
    plt.ylabel('Temperature (K)')
    plt.subplot(513)
    plt.scatter(data_list[0,0:-1:60],data_list[3,0:-1:60], s=4)
    plt.ylabel('Volume (LPM)')
    plt.subplot(514)
    plt.scatter(data_list[0,0:-1:60],data_list[4,0:-1:60], s=4)
    plt.ylabel('Mass (MPM)')
    plt.subplot(515)
    plt.scatter(data_list[0,0:-1:60],data_list[5,0:-1:60], s=4)
    plt.ylabel('Mass (kgpm)')
    plt.xlabel('Time (hours)')


#    plt.scatter(data_list[0,:],data_list[2,:], s=4)
#    ax=plt.gca()
#    ax.set_yscale('log')

#    ax = plt.gca()
#    ax.set_ylabel('Temperature ($ ^\\alpha $C)')
#    fig.savefig('all_vs_time.png', bbox_inches='tight')
#    plt.tight_layout()
    plt.tight_layout()
    fig.savefig('all_vs_time.png')
#    plt.show()

#def integrate_num_flux(seq):
    

def plot_flow_magnet_on(seq):
    data_list=load_file(1);
#    plt.yscale('log')
    fig=plt.figure(figsize=(10, 8), dpi=800,)
    start=2170 
    end=239030
    plt.subplot(411)
    plt.scatter(data_list[0,start:end],data_list[1,start:end],s=4)
    plt.title('Magnet Thermal Test, magnet on only')
    plt.ylabel('Pressure (PSI)')
    plt.subplot(412)
    plt.scatter(data_list[0,start:end],data_list[2,start:end],s=4)
    plt.ylabel('Temperature (K)')
    plt.subplot(413)
    plt.scatter(data_list[0,start:end],data_list[3,start:end],s=4)
    plt.ylabel('Volume (LPM)')
    plt.subplot(414)
    plt.scatter(data_list[0,start:end],data_list[4,start:end],s=4)
    plt.ylabel('Mass (MPM)')
    plt.xlabel('Time (hours)')


#    plt.scatter(data_list[0,:],data_list[2,:], s=4)
#    ax=plt.gca()
#    ax.set_yscale('log')

#    ax = plt.gca()
#    ax.set_ylabel('Temperature ($ ^\\alpha $C)')
#    fig.savefig('all_vs_time.png', bbox_inches='tight')
#    plt.tight_layout()
    plt.tight_layout()
    fig.savefig('magnet_on_vs_time.png')
#    plt.show()

#def integrate_num_flux(seq):

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
#    sum=0
#    i=0
#    while i<len(data_list[0,:]):
#        if i!=(len(data_list[0,:])-1):
#           sum+=(data_list[0,i+1]-data_list[0,i])
#           print (data_list[0,i+1]-data_list[0,i])
#        i+=1 
#    print (sum)
#    print(sum/len(data_list[0,:]))
    return data_list

def spline_levels_near(seq):
# get the data
    data_list=load_levels_near(1)
# convert back to minutes
#    time_mins= numpy.
#    data_list=numpy.array([time_mins,data_list[1,:]])
# spline section
    #interpolate the data using spline of order 3 (cubic by default)
    tck, fp, ier, msg= interpolate.splrep(data_list[0,:], data_list[1,:], s=0, full_output=True)
    #new array for x-axis values to evaluate the calculated spline at
    time_spline=numpy.arange(0, data_list[0,-1], 0.1)
    #get the yvalues from the found spline at the xvalues from the new array above. 
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

    
def load_levels_far(seq):
    datafar=numpy.genfromtxt('lvlSensorFar.csv', dtype=float, delimiter=',', names=True)
    timefar=datafar['time']
    lvlfar=datafar['NonStackSideLevelcm']
# convert unix time to seconds after tests started
    b=numpy.empty(len(timefar))
    b.fill(timefar[0])
    timefar=timefar-b    
# convert to minutes after test started
    timefar= numpy.true_divide(timefar,60.0)
    data_list=numpy.array([timefar,lvlfar])
    return data_list

def plot_levels_all(seq):
    data_listnear=load_levels_near(1);
    data_listfar=load_levels_far(1);
    fig=plt.figure(figsize=(10, 8), dpi=800)
#    plt.scatter(data_list[0,:],data_list[1,:],s=4)
    plt.scatter(data_listnear[0,:], data_listnear[1,:], c='b', marker='s', label='Near Sensor')
    plt.scatter(data_listfar[0,:], data_listfar[1,:], c='r', marker='s', label='Far Sensor')
#    ax1 = fig.add_subplot(111)
#    ax1.scatter(data_list[0,:], data_list[1,:], s=4, c='b', marker="s", label='near sensor')
#    ax1.scatter(data_list[0,:],data_list[2,:], s=4, c='r', marker="o", label='far sensor')

    plt.legend(loc='upper left');
    plt.title('Magnet Thermal Test, near and far sensor levels')
    plt.ylabel('Level (cm)')
    plt.xlabel('Time (minutes)')
    fig.savefig('magnet_lvls_vs_time_mins.png')

def plot_far_vs_near(seq):
    data_listnear=load_levels_near(1); # both zeroed to the same time
    data_listfar=load_levels_far(1); # both zeroed to the same time
    #only plot those values that are the same times
    near_spliced=[]
    j=0
    while j<len(data_listfar[1,:]):
       i=0
       while i<len(data_listnear[1,:]):
          if data_listnear[0,i]==data_listfar[0,j]: #match up the times
             near_spliced.append(data_listnear[1,i])
          i+=1
       j+=1
#    near_far=numpy.array([near_spliced,data_listfar[1,:]])
    print (len(near_spliced))
# find section also
    k=0
    lvlnear=[]
    lvlfar=[]
    while k< len(data_listfar[0,:]):
       if data_listfar[0,k]>2000 and data_listfar[0,k] <3900:
          lvlfar.append(data_listfar[1,k])
          lvlnear.append(near_spliced[k])
       k+=1
    fig=plt.figure(figsize=(10, 8), dpi=800)
    plt.scatter(near_spliced, data_listfar[1,:], c='b', s=3)
    plt.title('Magnet Thermal Test, near and far sensor levels')
    plt.ylabel('Far (cm)')
    plt.xlabel('Near(cm)')
    fig.savefig('magnet_lvls_far_vs_near_all.png')
#clear it
    plt.clf()
    plt.scatter(lvlnear, lvlfar, c='b', s=3)
    plt.title('Magnet Thermal Test, near and far sensor levels, section only')
    plt.ylabel('Far (cm)')
    plt.xlabel('Near(cm)')
    fig.savefig('magnet_lvls_far_vs_near_section.png')
# fit?
# try linear regression
    #do this tonight. 

        
def plot_spline_der_vs_near(seq):
    lvlnear=load_levels_near(1)
    # spline it and get derivative
    tck, fp, ier, msg= interpolate.splrep(lvlnear[0,:], lvlnear[1,:], s=0, full_output=True)
    #new array for x-axis values to evaluate the calculated spline at
#    time_spline=numpy.arange(0, lvlnear[0,-1], 0.1)
    time_spline=lvlnear[0,:]
    lvl_spline=interpolate.splev(time_spline, tck, der=0)
    lvl_der_spline=interpolate.splev(time_spline, tck, der=1)    
    # find the spline values only for times in which we have data points
#    der_spliced=[]
#    j=0
#    while j<len(lvlnear[1,:]):
#       i=0
#       while i<len(time_spline):
#          if lvlnear[0,j]==time_spline[i]: #match up the times
#             der_spliced.append(lvl_der_spline[i])
#          i+=1
#       j+=1
#lvlnear is bigger by 1 than der_spliced
    # plot the spline vs the lvl
    # find section 
    k=0
    nearsection=[]
    nearsectionder=[]
    while k< len(lvlnear[0,:]):
       if lvlnear[0,k]>2000 and lvlnear[0,k] <3900:
          nearsection.append(lvlnear[1,k])
          nearsectionder.append(lvl_der_spline[k])
       k+=1
    fig=plt.figure(figsize=(10, 8), dpi=800)
    plt.scatter(lvlnear[1,:], lvl_der_spline, c='b', s=3)
    plt.title('Magnet Thermal Test, near derivative and near sensor level')
    plt.ylabel('Near derivative (cm per min)')
    plt.xlabel('Near(cm)')
    fig.savefig('magnet_near_der_vs_near.png')
    plt.clf()
    plt.scatter(nearsection, nearsectionder, c='b', s=3)
    plt.title('Magnet Thermal Test, near derivative and near sensor level, Section only')
    plt.ylabel('Near derivative (cm per min)')
    plt.xlabel('Near(cm)')
    fig.savefig('magnet_near_der_vs_near_section.png')

def plot_flow_vs_near(seq):
    flows, lvlnear=load_flows_zero_to_near_lvl_sensor(1)
    sum=0
    m=0
    while m<len(lvlnear[0,:]):
       if lvlnear[0,m]<flows[0,0]:
           sum+=1
       m+=1
    print (sum)
    #match up times? x values are lvlnear[1,:] and do not change. so find the elements i of flows[3,i] 
    flow_spliced=[]
    lvl_spliced=[]
    j=0
    while j<len(lvlnear[0,:]):
       i=0
       while i<len(flows[0,:]):
          if lvlnear[0,j]>2000 and lvlnear[0,j]<3900:
             if lvlnear[0,j]>flows[0,i] and lvlnear[0,j]<flows[0,i+1]: #average the values nearest the right time? by the minute?
                average=(flows[3,i]+flows[3,i+1])/2.0
                flow_spliced.append(average) # the value nearest the time maybe?
                lvl_spliced.append(lvlnear[1,j])
          i+=1
       j+=1
    print(len(flow_spliced))
    print(len(lvl_spliced))
    fig=plt.figure(figsize=(10, 8), dpi=800)
    plt.scatter(lvl_spliced, flow_spliced, c='b', s=3)
    plt.title('Magnet Thermal Test, Flow and near sensor level')
    plt.ylabel('Flow (liters per min)')
    plt.xlabel('Near(cm)')
#plt.savefig('destination_path.eps', format='eps', dpi=1000)
    fig.savefig('magnet_flow_vs_near_section.eps')
#    plt.clf()
#    plt.scatter(lvl_spliced, flow_spliced, c='b', s=3)
#    plt.title('Magnet Thermal Test, Flow and near sensor level, section only')
#    plt.ylabel('Flow (liters per min)')
#    plt.xlabel('Near(cm)')
#    fig.savefig('magnet_flow_vs_near_section.png')

def plot_flow_over_der_vs_near(seq):
    # only wants array elements that are in all three of, flows, lvlnear, and der.
    flows, lvlnear=load_flows_zero_to_near_lvl_sensor(1)
# match up the flows to the lvlnear times
    flow_spliced=[]
    mass_flow_spliced=[]
    lvl_spliced=[]
    lvl_spliced_time=[]
    j=0
    while j<len(lvlnear[0,:]):
       i=0
       while i<len(flows[0,:]):
#          if lvlnear[0,j]>2000 and lvlnear[0,j]<3900:
          if lvlnear[0,j]>flows[0,i] and lvlnear[0,j]<flows[0,i+1]: #average the values nearest the right time? by the minute?
             average=(flows[3,i]+flows[3,i+1])/2.0
             average_mass=(flows[5,i]+flows[5,i+1])/2.0
             flow_spliced.append(average) # the value nearest the time maybe?
             mass_flow_spliced.append(average_mass) # the value nearest the time maybe?
             lvl_spliced.append(lvlnear[1,j])
             lvl_spliced_time.append(lvlnear[0,j])
          i+=1
       j+=1
    print(len(flow_spliced))
    print(len(lvl_spliced))
   # get derivative
    tck, fp, ier, msg= interpolate.splrep(lvl_spliced_time, lvl_spliced, s=0, full_output=True)
    #new array for x-axis values to evaluate the calculated spline at
#    time_spline=numpy.arange(0, lvlnear[0,-1], 0.1)
    time_spline=lvl_spliced_time
    lvl_spline=interpolate.splev(time_spline, tck, der=0)
    lvl_der_spline=interpolate.splev(time_spline, tck, der=1)    
    ratio=numpy.true_divide(flow_spliced,lvl_der_spline)
    ratio_mass=numpy.true_divide(mass_flow_spliced,lvl_der_spline) 
 # check some of the values?
    p=0
    ratio_spliced=[]
    ratio_mass_spliced=[]
    lvl_section=[]
    while p<len(ratio):
       if ratio[p]<800 and ratio[p]>-800:
          lvl_section.append(lvl_spliced[p])
          ratio_spliced.append(ratio[p])
          ratio_mass_spliced.append(ratio_mass[p])
       p+=1
    fig=plt.figure(figsize=(10, 8), dpi=800)
    plt.scatter(lvl_section, ratio_spliced, c='b', s=3)
    plt.title('Magnet Thermal Test, Flow over derivative and near sensor level, section only')
    plt.ylabel('Flow over derivative (liters per cm)')
    plt.xlabel('Near(cm)')
#plot the sections 
    fig.savefig('magnet_flow_over_der_vs_near_section.eps')
    plt.clf()
    plt.scatter(lvl_section, ratio_mass_spliced, c='b', s=3)
    plt.title('Magnet Thermal Test, Mass Flow over derivative and near sensor level, section only')
    plt.ylabel('Mass Flow over derivative (g per cm)')
    plt.xlabel('Near(cm)')
    fig.savefig('magnet_mass_flow_over_der_vs_near_section.eps')
# plot the entire thing
    plt.clf()
    plt.scatter(lvl_spliced, ratio, c='b', s=3)
    plt.title('Magnet Thermal Test, Mass Flow over derivative and near sensor level')
    plt.ylabel('Mass Flow over derivative (g per cm)')
    plt.xlabel('Near(cm)')
    fig.savefig('magnet_flow_over_der_vs_near.eps')
    # mass
    plt.clf()
    plt.scatter(lvl_spliced, ratio_mass, c='b', s=3)
    plt.title('Magnet Thermal Test, Mass Flow over derivative and near sensor level')
    plt.ylabel('Mass Flow over derivative (g per cm)')
    plt.xlabel('Near(cm)')
    fig.savefig('magnet_mass_flow_over_der_vs_near.eps')

   
def plot_flow_and_levels(seq):
# load lvlnear data
    datanear=numpy.genfromtxt('lvlSensorNear.csv', dtype=float, delimiter=',', names=True)
    timenear = datanear['time']
    lvlnear=datanear['StackSideLevelcm']
#load lvlfar data
    datafar=numpy.genfromtxt('lvlSensorFar.csv', dtype=float, delimiter=',', names=True)
    timefar=datafar['time']
    lvlfar=datafar['NonStackSideLevelcm']
#load flow data
    data=numpy.genfromtxt('WhisperData_Cleaned.csv', dtype=float, delimiter=',', names=True)
    timeflow = data['time_s']
    pressure = data['pressure_PSI']
    temp = data['temperatureC']
    flow = data['volume_flow_LPM']
    mass = data['mass_flow']
#load event times
    datanotes=numpy.genfromtxt('notable_LHe_events.csv', delimiter=',', names=True)
    xposition=datanotes['time_event']
    labels=datanotes['note']
# convert unix time to seconds after tests started for all arrays of time, each array of time has different lengths since data was poorly recorded.
   # make array of same length as each time array
    a=numpy.empty(len(timeflow))
    b=numpy.empty(len(timenear))
    c=numpy.empty(len(timefar))
    d=numpy.empty(len(xposition))
   # fill each array of time with entry of earliest unix time (timenear[0])
    a.fill(timenear[0])
    b.fill(timenear[0])
    c.fill(timenear[0])
    d.fill(timenear[0])
  # subtract each start time array (a,b,c) from the corresponding time array.
    timeflow=timeflow-a
    timenear=timenear-b
    timefar=timefar-c
    xposition=xposition-d
# convert to hours after test started
    timeflow= numpy.true_divide(timeflow,3600.0)
    timenear= numpy.true_divide(timenear,3600.0)
    timefar= numpy.true_divide(timefar,3600.0)
    xposition= numpy.true_divide(xposition,3600.0)
#calculate the hour marks for notable events for vertical lines
#    time1=1510279237
#    print xposition
#    xposition = [data_list[0,time0], data_list[0,time2], data_list[0,time1], data_list[0,-1]] # these will be values in hours.
#    labels=['Ramp start, FULL', 'start ramping ','Ramp end, FULL','end of data']

# one big numpy array for the flows
    data_list=numpy.array([timeflow,pressure,temp,flow,mass])
# array for lvl near
    data_near=numpy.array([timenear,lvlnear])
# array for lvl far
    data_far=numpy.array([timefar,lvlfar])

#integrate the volume flow to get total volume flowed out. 
    #need to think if this is for entire time or part of time
#mark on graphs when LHe was added ~1 hour or so to transfer, "topping off" procedure
 
    #first, mark on the graph where those LHE times are....
    fig=plt.figure(figsize=(10, 8), dpi=800,)

    ax_flow=plt.subplot(311)
    ax_flow.scatter(data_list[0,:],data_list[3,:],s=4)
#put the lines for the test details here
# need to put these lines (but not the labels) on each subplot. Need to resize the text and add more lines for other notable things during data taking. 
    i=0
    for xc in xposition:
      ax_flow.axvline(x=xc, color='k', linestyle='--') # this is for the lines to mark at what time notable things in the test occured.
#      ax_flow.text(xc+1, 15,labels[i], rotation=90, fontsize=6)
      i+=1
    ax_flow.set_xlim([0,300])

    ax_flow.set_ylabel('Volume (LPM)')

    ax_near=plt.subplot(312)
    ax_near.scatter(data_near[0,:],data_near[1,:],s=4)
    ax_near.set_ylabel('Level Near(cm)')
    ax_near.set_xlim([0,300])

    ax_far=plt.subplot(313)
    ax_far.scatter(data_far[0,:],data_far[1,:],s=4)
    ax_far.set_ylabel('Level Far(cm)')
    ax_far.set_xlim([0,300])

    ax_far.set_xlabel('Time (hours)')
#    plt.tight_layout()
    fig.savefig('magnet_flowmeterquantities_vs_time.png')

def load_flows_zero_to_near_lvl_sensor(seq):
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

def plot_lvl_near_sensor_vs_time_section(seq):
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
        if flows[0,k] >5138 and flows[0,k]<6826.13 and flows[3,k]>20 and flows[3,k]<28:
           flows_spliceA.append(flows[3,k])
           time_spliceA.append(flows[0,k])
           #how to splice simultaneously the lvl sensor? its a smaller array so just go through the elements?
        if flows[0,k]>6900 and flows[0,k]<8382.8:
           flows_spliceB.append(flows[3,k])
           time_spliceB.append(flows[0,k])
        if flows[0,k]>10000 and flows[0,k]<16500:
           flows_spliceC.append(flows[3,k])
           time_spliceC.append(flows[0,k])
        k+=1
###splice levels into sections###
    j=0
    while j<len(lvlnear[0,:]):
        if lvlnear[0,j] >time_spliceA[0] and lvlnear[0,j]<time_spliceA[-1]:
           lvl_spliceA.append(lvlnear[1,j])
           lvl_time_spliceA.append(lvlnear[0,j])
        if lvlnear[0,j] >time_spliceB[0] and lvlnear[0,j]<time_spliceB[-1]:
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
    fig=plt.figure(figsize=(10, 8), dpi=800)
    if seq==1:
       fig.subplots_adjust(hspace=0)
       #need subplots, subplots(sharex=True)
       ax_flow=plt.subplot(211)
       plt.setp(ax_flow.get_xticklabels(), visible=False)
       ax_flow.scatter(flows[0,:],flows[3,:],s=3)
       # add in the section lines
       xposition=[5138,6826.13,6900,8382.8,10000,16500]
       labels=['section A', 'section B', 'section C']
       label_pos=[5300,7000,12000]
       i=0
       for xc in xposition:
          ax_flow.axvline(x=xc, color='k', linestyle='--') # this is for the lines to mark at what time notable things in the test occured.
       for xl in labels:
          ax_flow.text(label_pos[i], 40, labels[i], rotation=0, fontsize=8)
          i+=1
       ax_flow.set_ylabel('Flow (liters per minute)')
       ax_flow.set_xlim([0,18000])
       ax_near=plt.subplot(212,sharex=ax_flow)
       ax_near.scatter(lvlnear[0,:],lvlnear[1,:],s=3)
       for xi in xposition:
          ax_near.axvline(x=xi, color='k', linestyle='--') # this is for the lines to mark at what time notable things in the test occured.
       ax_near.set_ylabel('Level Sensor (cm)')
       ax_near.set_xlim([0,18000])
       ax_near.set_xlabel('Time (mins)')
       fig.savefig('ALL_flow_and_lvl_vs_time.eps')
       plt.clf()
#################################################SECTIONS#############################################################
###plot sectionA###
    # plot flow
    if seq==2:
       plt.scatter(time_spliceA, flows_spliceA, c='b', s=3)
       plt.title('Flow vs Time, section A')
       plt.ylabel('Flow (liters per min)')
       plt.xlabel('Time (mins)')
       fig.savefig('sectionA_flow_vs_time.eps')
       plt.clf()
       # plot level data and spline
   #    plt.scatter(lvl_time_spliceA, splineA, c='r', s=4, marker='s', label='Spline')
       plt.scatter(splineA_times_more, splineA_more, c='r', s=4, marker='s', label='Spline')
       plt.scatter(lvl_time_spliceA, lvl_spliceA, c='b', s=3, label='Data')
       plt.legend(loc='upper right')
       plt.title('Level vs Time, section A')
       plt.ylabel('Level (cm)')
       plt.xlabel('Time (mins)')
       fig.savefig('sectionA_lvl_and_spline_vs_time.eps')
   #    fig.savefig('sectionA_lvl_vs_time.eps')
       plt.clf()
      #plot flow and level with subplots?
   #https://matplotlib.org/devdocs/gallery/subplots_axes_and_figures/shared_axis_demo.html#sphx-glr-gallery-subplots-axes-and-figures-shared-axis-demo-py
   #    ax_spline=plt.subplot(211,sharex
    if seq==1:
# make subplots instead 
       fig.subplots_adjust(hspace=0)
       ax_Aflow=plt.subplot(211)
       ax_Aflow.set_title('Section A')
       plt.setp(ax_Aflow.get_xticklabels(), visible=False)              
       ax_Aflow.scatter(time_spliceA, flows_spliceA, c='b', s=3)
       ax_Aflow.set_ylabel('Flow (liters per min)')
       ax_Aspline=plt.subplot(212)
       ax_Aspline.scatter(splineA_times_more, splineA_der_more, c='r', s=4, marker='s', label='Spline')
       ax_Aspline.scatter(lvl_time_spliceA, splineA_der, c='b', s=4, marker='s', label='Spline at data points')
#       ax_Aspline.scatter(lvl_time_spliceA, lvl_spliceA, c='b', s=3, label='Data')
       ax_Aspline.legend(loc='upper right')
       ax_Aspline.set_ylabel('Level derivative (cm per min)')
       ax_Aspline.set_xlabel('Time (mins)')
       fig.savefig('sectionA_lvl_der_spline_vs_time_subplots_both.eps')
#       fig.savefig('sectionA_lvl_der_spline_vs_time_subplots_more.eps')
       plt.clf()
    if seq==3:
    #plot flows and derivative at all flow points
       fig.subplots_adjust(hspace=0)
       ax_Aflow=plt.subplot(211)
       ax_Aflow.set_title('Section A')
       plt.setp(ax_Aflow.get_xticklabels(), visible=False)              
       ax_Aflow.scatter(time_spliceA, flows_spliceA, c='b', s=3)
       ax_Aflow.set_ylabel('Flow (liters per min)')
       ax_Aspline=plt.subplot(212)
       ax_Aspline.scatter(time_spliceA, splineA_der_flowtimes, c='r', s=4, marker='s', label='Spline')
#       ax_Aspline.scatter(lvl_time_spliceA, splineA_der, c='r', s=4, marker='s', label='Spline')
#       ax_Aspline.scatter(lvl_time_spliceA, lvl_spliceA, c='b', s=3, label='Data')
 #      ax_Aspline.legend(loc='upper right')
       ax_Aspline.set_ylabel('Level derivative (cm per min)')
       ax_Aspline.set_xlabel('Time (mins)')
       fig.savefig('sectionA_lvl_der_spline_vs_time_subplots_flowtimes.eps')
       plt.clf()
#    if seq==4:
       #make plot of ratios, which involves calculating the flows at the lvltimes?

###plot sectionB###
    if seq==2:
       plt.scatter(time_spliceB, flows_spliceB, c='b', s=3)
       plt.title('Flow vs Time, section B')
       plt.ylabel('Flow (liters per min)')
       plt.xlabel('Time (mins)')
       fig.savefig('sectionB_flow_vs_time.eps')
       plt.clf()
#    plt.scatter(lvl_time_spliceB, lvl_spliceB, c='b', s=3)
       plt.scatter(splineB_times_more, splineB_more, c='r', s=4, marker='s', label='Spline')
       plt.scatter(lvl_time_spliceB, lvl_spliceB, c='b', s=3, label='Data')
       plt.legend(loc='upper right')
       plt.title('Level vs Time, section B')
       plt.ylabel('Level (cm)')
       plt.xlabel('Time (mins)')
#       fig.savefig('sectionB_lvl_and_spline_vs_time.eps')
#    fig.savefig('sectionB_lvl_vs_time.eps')
       plt.clf()
    if seq==1:
# make subplots instead 
       fig.subplots_adjust(hspace=0)
       ax_Bflow=plt.subplot(211)
       ax_Bflow.set_title('Section B')
       plt.setp(ax_Bflow.get_xticklabels(), visible=False)              
       ax_Bflow.scatter(time_spliceB, flows_spliceB, c='b', s=3)
       ax_Bflow.set_ylabel('Flow (liters per min)')
       ax_Bspline=plt.subplot(212)
       ax_Bspline.scatter(splineB_times_more, splineB_der_more, c='r', s=4, marker='s', label='Spline')
       ax_Bspline.scatter(lvl_time_spliceB, splineB_der, c='b', s=4, marker='s', label='Spline at data points')
#       ax_Bspline.scatter(lvl_time_spliceB, lvl_spliceB, c='b', s=3, label='Data')
       ax_Bspline.legend(loc='upper right')
       ax_Bspline.set_ylabel('Level derivative (cm per min)')
       ax_Bspline.set_xlabel('Time (mins)')
       fig.savefig('sectionB_lvl_der_spline_vs_time_subplots_both.eps')
#       fig.savefig('sectionB_lvl_der_spline_vs_time_subplots_more.eps')
       plt.clf()
    if seq==3:
# make subplots instead 
       fig.subplots_adjust(hspace=0)
       ax_Bflow=plt.subplot(211)
       ax_Bflow.set_title('Section B')
       plt.setp(ax_Bflow.get_xticklabels(), visible=False)              
       ax_Bflow.scatter(time_spliceB, flows_spliceB, c='b', s=3)
       ax_Bflow.set_ylabel('Flow (liters per min)')
       ax_Bspline=plt.subplot(212)
       ax_Bspline.scatter(time_spliceB, splineB_der_flowtimes, c='r', s=4, marker='s', label='Spline')
#       ax_Bspline.scatter(lvl_time_spliceB, splineB_der, c='r', s=4, marker='s', label='Spline')
#       ax_Bspline.scatter(lvl_time_spliceB, lvl_spliceB, c='b', s=3, label='Data')
#       ax_Bspline.legend(loc='upper right')
       ax_Bspline.set_ylabel('Level derivative (cm per min)')
       ax_Bspline.set_xlabel('Time (mins)')
       fig.savefig('sectionB_lvl_der_spline_vs_time_subplots_flowtimes.eps')
       plt.clf()

###plot sectionC###
    if seq==2:
       plt.scatter(time_spliceC, flows_spliceC, c='b', s=3)
       plt.title('Flow vs Time, section C')
       plt.ylabel('Flow (liters per min)')
       plt.xlabel('Time (mins)')
       fig.savefig('sectionC_flow_vs_time.eps')
       plt.clf()
#       plt.scatter(lvl_time_spliceC, lvl_spliceC, c='b', s=3)
       plt.scatter(splineC_times_more, splineC_more, c='r', s=4, marker='s', label='Spline')
       plt.scatter(lvl_time_spliceC, lvl_spliceC, c='b', s=3, label='Data')
       plt.legend(loc='upper right')
       plt.title('Level vs Time, section C')
       plt.ylabel('Level (cm)')
       plt.xlabel('Time (mins)')
       fig.savefig('sectionC_lvl_and_spline_vs_time.eps')
#       fig.savefig('sectionC_lvl_vs_time.eps')
    if seq==1:
# make subplots instead 
       fig.subplots_adjust(hspace=0)
       ax_Cflow=plt.subplot(211)
       ax_Cflow.set_title('Section C')
       plt.setp(ax_Cflow.get_xticklabels(), visible=False)              
       ax_Cflow.scatter(time_spliceC, flows_spliceC, c='b', s=3)
       ax_Cflow.set_ylabel('Flow (liters per min)')
       ax_Cspline=plt.subplot(212)
       ax_Cspline.scatter(splineC_times_more, splineC_der_more, c='r', s=4, marker='s', label='Spline')
       ax_Cspline.scatter(lvl_time_spliceC, splineC_der, c='b', s=4, marker='s', label='Spline at data points')
#       ax_Cspline.scatter(lvl_time_spliceC, lvl_spliceC, c='b', s=3, label='Data')
       ax_Cspline.legend(loc='upper right')
       ax_Cspline.set_ylabel('Level derivative (cm per min)')
       ax_Cspline.set_xlabel('Time (mins)')
       fig.savefig('sectionC_lvl_der_spline_vs_time_subplots_both.eps')
#       fig.savefig('sectionC_lvl_der_spline_vs_time_subplots_more.eps')
#       fig.savefig('sectionC_lvl_der_spline_vs_time_subplots.eps')
#       plt.clf()
    if seq==3:
# make subplots instead 
       fig.subplots_adjust(hspace=0)
       ax_Cflow=plt.subplot(211)
       ax_Cflow.set_title('Section C')
       plt.setp(ax_Cflow.get_xticklabels(), visible=False)              
       ax_Cflow.scatter(time_spliceC, flows_spliceC, c='b', s=3)
       ax_Cflow.set_ylabel('Flow (liters per min)')
       ax_Cspline=plt.subplot(212)
       ax_Cspline.scatter(time_spliceC, splineC_der_flowtimes, c='r', s=4, marker='s', label='Spline')
#       ax_Cspline.scatter(lvl_time_spliceC, splineC_der_more, c='r', s=4, marker='s', label='Spline')
#       ax_Cspline.scatter(lvl_time_spliceC, lvl_spliceC, c='b', s=3, label='Data')
#       ax_Cspline.legend(loc='upper right')
       ax_Cspline.set_ylabel('Level derivative (cm per min)')
       ax_Cspline.set_xlabel('Time (mins)')
       fig.savefig('sectionC_lvl_der_spline_vs_time_subplots_flowtimes.eps')
###################PLOT RATIOS######################
# think about moving averages like this 
#def moving_average(a, n=3) :
#    ret = np.cumsum(a, dtype=float)
#    ret[n:] = ret[n:] - ret[:-n]
#    return ret[n - 1:] / n

##### calculate ratio in more than one scenario####
##moving average


