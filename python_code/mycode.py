filename="Ratio_vs_levels_binwidth_subplots.png"
    fig.savefig(filename)
    plt.clf()
    cryo_data=numpy.genfromtxt('vol_vs_lvl_cryo.csv', dtype=float, delimiter=',', names=True)
    lvl_cryo = cryo_data['lvl']
    vol_cryo = cryo_data['vol']
    
    
    #       plt.scatter(lvl_cryo,vol_cryo,s=3, c='black')
    x=numpy.arange(lvl_cryo[2],100,0.01)
        #convert fits to liters (multiply by 10000 for cm^2 then divide by 1000 for the conversion to liters) and add in the value of cryo at around 30cm (since data stops there?)
        fit40=(pvals40(50)*10.0)*x+(vol_cryo[2]-(pvals40(50)*10.0*x[0]))
        fit20=(pvals20(50)*10.0)*x+(vol_cryo[2]-(pvals20(50)*10.0*x[0]))
        fit10=(pvals10(50)*10.0)*x+(vol_cryo[2]-(pvals10(50)*10.0*x[0]))
        fit5=(pvals5(50)*10.0)*x+(vol_cryo[2]-(pvals5(50)*10.0*x[0]))
        plt.plot(x,fit40,'-.',x,fit20,':',x,fit10,x,fit5,'--',lvl_cryo,vol_cryo,'o')
        plt.xlim(30,100)
        plt.ylim(0,300)
        
        nvol=[]
        nlvl=[]
        j = 0
        while j < len(lvl_cryo):
            if lvl_cryo[j]<=100 and lvl_cryo[j]>=30:
                nlvl.append(lvl_cryo[j])
                nvol.append(vol_cryo[j])
                print(nlvl[-1])
            j+=1
                            
                            
         
                                    
                                    
       x2 =numpy.arange(nlvl[0],nlvl[-1],0.01)
       y2=numpy.polyfit(nlvl,nvol,1)
       y2[1]=y2[1]+nvol[0]
       z2=numpy.polyval(y2,x2)
       plt.plot(x2,z2,c='red',dashes=[6,2])
       plt.xlim(30,100)
       plt.ylim(0,300)
                                    
        slope  = format(pvals40(50)*10.0,'4.3f')
        slope2 = format(pvals20(50)*10.0,'4.3f')
        slope3 = format(pvals10(50)*10.0,'4.3f')
        slope4 = format(pvals5(50)*10.0,'4.3f')
        slope5 = format(y2[0],'4.3f')
                                    
                                    
         plt.legend(("m="+slope+"  ∆L=1.5375","m="+slope2+"  ∆L=3.075","m="+slope3+"  ∆L=6.15","m="+slope4+"  ∆L=12.3",'Data',"m="+slope5+"  Fit data"),loc='upper left',prop={'size': 18})
         plt.title('Geometric calibration using gaseous flow',fontsize=20)
         plt.ylabel('Volume (liters)',fontsize=16)
         plt.xlabel('Sensor Level (cm)',fontsize=16)
                                    
         plt.tick_params(labelsize=16)
                                    
         filename="Volume_vs_levels_cryo_doc.png"
         fig.savefig(filename)
