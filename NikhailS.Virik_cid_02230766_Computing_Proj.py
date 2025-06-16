# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 16:36:40 2022

@author: Nikhail
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
lest=list()
lest2=list() #lists to append oberved frequency later on


l_data=np.loadtxt(r'C:\Users\Nikhail\Downloads\Data_proj\Datap.2\Halpha_spectral_data.csv', skiprows=3, delimiter=',', unpack=True) #l_data holds the h-alpha data
dist_data=np.loadtxt(r'C:\Users\Nikhail\Downloads\Data_proj\Datap.2\Distance_Mpc.txt', delimiter='\t', skiprows=2) #dist_data holds the distance data
row_1=l_data[:,0]
all_id=row_1[::2] #parses out all the observation numbers the the fourth row of h-alpha (obsrv nums) and arrays in all_id
all_id=all_id.flatten() #np.flatten() code was found at https://numpy.org/, Author: Charis et al, "numpy.ndarray.flatten", [Online]. Available: https://numpy.org/doc/stable/reference/generated/numpy.ndarray.flatten.html, Accessed on: Nov. 11, 2022
d_id=dist_data[:,:1] #parses out the observation numbers from the distance data and arrays it into d_id
d_id=d_id.flatten() #np.flatten() code was found at https://numpy.org/, Author: Charis et al, "numpy.ndarray.flatten", [Online]. Available: https://numpy.org/doc/stable/reference/generated/numpy.ndarray.flatten.html, Accessed on: Nov. 11, 2022
dist=dist_data[:,1:2] #parses out all the distances and arrays them in array named dist
dist=dist.flatten()
dist=np.array(dist)
valid=dist_data[:,2:] #parses out the instrument responses and arrays them in array named valid
valid=valid.flatten() ##np.flatten() code was found at https://numpy.org/, Author: Charis et al, "numpy.ndarray.flatten", [Online]. Available: https://numpy.org/doc/stable/reference/generated/numpy.ndarray.flatten.html, Accessed on: Nov. 11, 2022
id_id=l_data[:,0]
id_1=id_id[::2] #creates an array of the observation numbers from h-alpha and removes the double entry. This is arrayed into id_id





for i in range(0,(len(l_data)),2):
    set_1=l_data[i:i+2] #iterates thorugh each set of observation data in h-alpha
    set_id=set_1[:,0][0] #iterates through the observation numbers in h-alpha. Obsrv nums are stored in array named set_id
    frq_1=set_1[:,1:][0] #interates through the frequncy data in h-alpha. Frequencies are stored in array called frq_1
    lum_1=set_1[:,1:][1] #interates through the spectra data in h-alpha. Spectral data is stored in array called lum_1
    plt.plot(frq_1,lum_1) #plots the curves of spectra against frequency
    
    def fit_func(frq_1,a,mu,sig,m,c):
        gauss=a*np.exp(-((frq_1)-mu)**2/(2*sig**2)) # formula for a gaussian curve; based on the frequency data, a is the amplitude of the gaussian, mu is the mean value of the gaussian, sig is the std deviation width
        line=m*(frq_1)+c # formula for the straight line, m is the gradient, c is the y-int
        return gauss + line #returns the gaussian and straight line plotted together
  
    
    grad_guess=((lum_1[-1]-lum_1[0])/(frq_1[-1]-frq_1[0])) # creates a guess of gradient of straight line by gradient formula on first and last (x,y) coordinates
                 
    if i==38: #correction for curve with label 38 with very low dip just before frq=4.5e14 such that maximum deviation from the gradient is not on the gaussian
        rang=np.max(lum_1)-np.min(lum_1)
    else:
        rang=np.max(lum_1)-grad_guess #a guess for the amplitude of the gaussian made on the deviation from the standard gradient
    ig=[rang,np.mean(frq_1),np.std(frq_1, ddof=1),grad_guess,lum_1[0]] #array of initial guesses. initial guess for mean of the gaussian was the mean of the frequency data, guess of width of gaussian was std dev in freq data, initial guess for y-int was the intensity at the smallest freq in the data
    try:
        po,po_cov=curve_fit(fit_func,frq_1,lum_1,ig)
        plt.plot(frq_1,fit_func(frq_1,po[0],po[1],po[2],po[3],po[4])) #corresponding returned values in the same order as our initial guesses
        plt.ylabel('Intensity(a.u)')
        plt.xlabel('frequency(Hz)')
        plt.legend([i], loc=1) #legend made for easily identifying curves
        plt.show()
        mean=po[1]
        # this code creates a straight line fit for the straight parts of the curves and a gaussian for the curved parts, iterating based on initial guesses until fit is acheieved
        # mean gauss freq was taken as observed freq. A for loop interated through each set of data and appened each observed freq into lest
        lest.append(mean)
        
        for val in lest:
            if val not in lest2:
                lest2.append(val)
        #a new list is created with every iteration. This code then takes unique values to create one list of all the observed freq      
        #Obsrvd freqs are now listed in lest2              
    except:
        continue
lest2=np.array(lest2) #lest2 is now an array of the obsrvd frqs

leste=list()
lestf=list()#lesftf contains all the data contianed within 0 instrument response, so they can be ignored in the data analysis

for p in range(0,len(id_1)):
    for q in range(0,len(d_id)):
        for r in range(0,len(valid)):
            if valid[r]==0: 
                x=d_id[r] #this code says that is the value in valid==0,find the obsrv num corresponding to this data set from d_id
                leste.append(x) #this code creates a list of all the observation numbers that correspond to bad instrument responses and lists them
                for bad in leste:
                    if bad not in lestf:
                        lestf.append(bad)
                    # this is a similar code to previous to condense all the lists into 1 list
                    else:
                        continue
            else:
                continue
             


alphlst=list()
bravlst=list()
for alpha in range(0,len(lestf)):
    for bravo in range(0,len(d_id)):
        if lestf[alpha]==d_id[bravo]:
            alphlst.append(bravo)
        #creates a list of the positions of the "bad" observation numbers in Distance_Mpc
alphlst=np.array(alphlst)      
newd_id=np.delete(d_id, alphlst) #np.delete code was found at https://geeksforgeeks.org, Aunthor: Mohit Gupta_OMG, aryen1713045(Ed.), kothavvsaakash(Ed.), vivekedule(Ed.), "numpy.delete() in Python" [Online]. Available: https://www.geeksforgeeks.org/numpy-delete-python/#:~:text=The%20numpy.,along%20with%20the%20mentioned%20axis.&text=Return%20%3A,object%20along%20a%20given%20axis., Accessed on: Nov. 17, 2022 
newdist=np.delete(dist, alphlst) #np.delete code was found at https://geeksforgeeks.org, Aunthor: Mohit Gupta_OMG, aryen1713045(Ed.), kothavvsaakash(Ed.), vivekedule(Ed.), "numpy.delete() in Python" [Online]. Available: https://www.geeksforgeeks.org/numpy-delete-python/#:~:text=The%20numpy.,along%20with%20the%20mentioned%20axis.&text=Return%20%3A,object%20along%20a%20given%20axis., Accessed on: Nov. 17, 2022
#this uses those positions to delete the "bad" observation numbers and corresponding distances from the d_id and dist arrays respectively, and arrays the new clean data
                
for charlie in range(0,len(lestf)):
    for delta in range(0,len(id_1)):
        if lestf[charlie]==id_1[delta]:
            bravlst.append(delta)
        #similarly, this creates a list of the positions of the "bad" observation numbers in the h-alpha data, remembering that the positions will be different for the two sets of data as the observation numbers are arranged differently in both
bravlst=np.array(bravlst)        
newid_1=np.delete(id_1,bravlst) #np.delete code was found at https://geeksforgeeks.org, Aunthor: Mohit Gupta_OMG, aryen1713045(Ed.), kothavvsaakash(Ed.), vivekedule(Ed.), "numpy.delete() in Python" [Online]. Available: https://www.geeksforgeeks.org/numpy-delete-python/#:~:text=The%20numpy.,along%20with%20the%20mentioned%20axis.&text=Return%20%3A,object%20along%20a%20given%20axis., Accessed on: Nov. 17, 2022
newlest2=np.delete(lest2, bravlst) #np.delete code was found at https://geeksforgeeks.org, Aunthor: Mohit Gupta_OMG, aryen1713045(Ed.), kothavvsaakash(Ed.), vivekedule(Ed.), "numpy.delete() in Python" [Online]. Available: https://www.geeksforgeeks.org/numpy-delete-python/#:~:text=The%20numpy.,along%20with%20the%20mentioned%20axis.&text=Return%20%3A,object%20along%20a%20given%20axis., Accessed on: Nov. 17, 2022 
# this code deletes the "bad" observation numbers and corresponding observed freqs and arrays the clean data

shift=list() #will hold redshift velocities
disti=list() #will hold distances(Mpc)
for o in range(0,len(newlest2)):
    
    y=(2.9979E8)/(newlest2[o])
    z=y/656.28E-9
    q=((2.9979E8)*((z*z)-1))/(1+(z*z))
    p=q/1000
    shift.append(p)
#this code iterates through the clean observed freqs and applies the formula to obtain recession speed from freq, and lists recession speeds into shift
shift=np.array(shift) #shift is then arrayed
lost=list()
for g in range(0,len(newid_1)):
    for h in range(0,len(newd_id)):
        if newid_1[g]==newd_id[h]:
            lost.append(newd_id[h])
            disti.append(newdist[h])
#this code iterates thorugh the clean obsrv nums in h-alpha and finds the same obsrv number in distance. The obsrv nums for ditance is then appended to a new list which is now in the same order as h-alpha called newd_id
#similarly, the corresponding clean distances that were originally found in position h in distance_mpc are put into position g which correlates them to the order of obsrv nums in h-alpha, essentially the data is moving from position g to h
newd_id=np.array(lost)  
disti=np.array(disti)       
#the newly arranged data is arrayed. disti holds the distance data in the correct order now
      
fit_shift,cov_shift=np.polyfit(disti,shift,1,cov=True) #polyfit for the data in shift, against distance data, returning covariance matrix wihtout pre-measured errors
pshift=np.poly1d(fit_shift) #allows us to return values from our graph 
plt.xlabel("Distance(Mpc)")
plt.ylabel("Redshift(km/s)")
plt.scatter(disti,shift)
plt.plot(disti,pshift(disti)) #this code plots the best fit-line on top of the scatter plots of data

plt.errorbar(disti,shift,xerr=0,yerr=np.std(np.sqrt(np.diag(cov_shift)),ddof=1),capsize=4, ls='none') #this code plots the error bars obtained from the curve fitting
#np.diag code was found at https://stackoverflow.com, Author: rikyborg, "Covariance matrix from np.polyfit() has negative diagonal?" [Online]. Available: https://stackoverflow.com/a/40104022, Accessed on: Nov. 21, 2022
plt.show()

print('H0=', fit_shift[0],"km/s/Mpc") #prints H0 
print('Uncertainty=', np.sqrt(cov_shift[0,0]),"km/s/Mpc")  #prints uncertainty in Ho           
print("Final Result: H0= %.3e +/- %.3e km/s/Mpc" %(fit_shift[0],np.sqrt(cov_shift[0,0])))     
plt.show()


        

