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
# convert to hours 
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

def plot_file(seq):
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
    
