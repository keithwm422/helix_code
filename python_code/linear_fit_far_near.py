#This code loads and fits the far level sensor and the near level sensor and computes residuals for the fits
import matplotlib.pyplot as plt
import numpy
from scipy import interpolate

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
# find linear section
    k=0
    lvlnear=[]
    lvlfar=[]
    while k< len(data_listfar[0,:]):
       if data_listfar[1,k]>20 and data_listfar[1,k] <95 and near_spliced[k]<110:
          lvlfar.append(data_listfar[1,k])
          lvlnear.append(near_spliced[k])
       k+=1
    fig=plt.figure(figsize=(10, 8), dpi=800)
    plt.scatter(near_spliced, data_listfar[1,:], c='b', s=3)
    plt.title('Magnet Thermal Test, near and far sensor levels')
    plt.ylabel('Far (cm)')
    plt.xlabel('Near(cm)')
    fig.savefig('lvls_far_vs_near_all.png')
#clear it
    plt.clf()
    plt.scatter(lvlnear, lvlfar, c='b', s=3)
    plt.title('Magnet Thermal Test, near and far sensor levels, section only')
    plt.ylabel('Far (cm)')
    plt.xlabel('Near(cm)')
    fig.savefig('lvls_far_vs_near_linear.png')
    plt.clf()
# fit?
# try linear regression
    #do this tonight. 
    p, res, rank, single,rcond=numpy.polyfit(lvlnear, lvlfar, 3, full=True)
    print (res)
    pxvals=numpy.linspace(20,100,100)
    pvals=numpy.poly1d(p)
    plt.scatter(lvlnear, lvlfar, color='b', s=20, label='Data')
    plt.scatter(pxvals, pvals(pxvals), color='g', s=10, label='polyfit')
    plt.title('Magnet Thermal Test, near and far sensor levels with fit')
    plt.legend(loc='upper left')
    plt.ylabel('Far (cm)')
    plt.xlabel('Near(cm)')
    fig.savefig('lvls_far_vs_near_with_fit.png')

