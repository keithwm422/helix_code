# Python code to plot the data from HELIX magnet test. Reads in .csv file downloaded form the HELIX wiki and plots versus time the quantities.  
# Author Keith McBride 10/16/18

import matplotlib.pyplot as plt
import numpy

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
    volume = data['volume_flow_LPM']
    mass = data['mass_flow']
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
    data_list=numpy.array([time,pressure,temp,volume,mass])
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
    data_list[0,:]= numpy.true_divide(data_list[0,:],0.016666666667)
    
    #integrate over time and report the value
    tot=numpy.trapz(data_list[4,:],x=data_list[0,:])
#    i=0
#    tot=0
#    while i<len(data_list[4,:]):
      
      
#      i+=1
    print('total mass flowed out is {}'.format(tot)) 

def integrate_volume_range(seq):
    data_list=load_file(1)
    fig=plt.figure(figsize=(10, 8), dpi=800,)
#    start = 2170 #found from grep -n "1510520760" WhisperCleaned_data.csv
#    end = 239030 # found from same as above but for time of 1510520760. also, because of 1st line in WhisperData_cleaned.csv being headers, subtract one from the line number. Then subtract one again for the accessing elements of an array (start at zero). 
    start = 239030 # found from same as above but for time of 1510520760. also, because of 1st line in WhisperData_cleaned.csv being headers, subtract one from the line number. Then subtract one again for the accessing elements of an array (start at zero). 
    end = len(data_list[0,:])-1 # found from same as above but for time of 1510520760. also, because of 1st line in WhisperData_cleaned.csv being headers, subtract one from the line number. Then subtract one again for the accessing elements of an array (start at zero). 

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
#    end = 239030 # found from same as above but for time of 1510520760. also, because of 1st line in WhisperData_cleaned.csv being headers, subtract one from the line number. Then subtract one again for the accessing elements of an array (start at zero). 
    start = 239030 # found from same as above but for time of 1510520760. also, because of 1st line in WhisperData_cleaned.csv being headers, subtract one from the line number. Then subtract one again for the accessing elements of an array (start at zero). 
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


def plot_file_all(seq):
    data_list=load_file(1);
#    plt.yscale('log')
    fig=plt.figure(figsize=(10, 8), dpi=800,)

    plt.subplot(411)
    plt.scatter(data_list[0,0:-1:60],data_list[1,0:-1:60], s=4)
    plt.title('Magnet Thermal Test')
    plt.ylabel('Pressure (PSI)')
    plt.subplot(412)
    plt.scatter(data_list[0,0:-1:60],data_list[2,0:-1:60], s=4)
    plt.ylabel('Temperature ($^\circ$C)')
    plt.subplot(413)
    plt.scatter(data_list[0,0:-1:60],data_list[3,0:-1:60], s=4)
    plt.ylabel('Volume (LPM)')
    plt.subplot(414)
    plt.scatter(data_list[0,0:-1:60],data_list[4,0:-1:60], s=4)
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
    fig.savefig('all_vs_time.png')
#    plt.show()

#def integrate_num_flux(seq):
    

def plot_file_magnet_on(seq):
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
    plt.ylabel('Temperature ($^\circ$C)')
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
    
