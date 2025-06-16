import numpy as np
import matplotlib.pyplot as plt
import STOM_higgs_tools as stom
n=0
total_chi=list()
total_p=list()
while n<=0:
    #for l in range(0,200):
                values = stom.generate_data() 
            
            #print(np.size(values))
            #We are taking an uncertainty of the y-values as the customary sqrt(n) where n is the number of values in each bin
            #and the bin width in the x-value errors
            #Create a numpy histogram so that we can categorise the photon masses into bins
            
            
    
                bin_heights, bin_edges = np.histogram(values, range=[100,160], bins=35) 
                bin_centre = bin_edges[1:] - 6/7 #Centre the bins
                
                
                '''
                Q1: Compare to figure 1 and how does the data visaully change with bin number?
                    Both graphs are extremely similar, they both appear to have a peak/bump between 120-130, with a little jump around the 
                135 mark (our graph shows the jumps a little earlier than the real data). They both appear to follow an exponential decay 
                curve with a small gaussian bump around the centre of the graph
                    The graph changes quite darastically with an increase in the number of bins. The data becomes smoother, and it is clear
                where the guassian bump appears, however after a certain amount of bins the points appear on top of each other making it harder
                to read the data (found this to be >35 bins)
                
                
                #The following Section is to calculate the lamda and chi values, however they take forever to output as it's not an efficient
                #algorithm, it outputs the following data:
                #Lamda = 30.0
                #A = 56560
                
                
                #For this section chose data where there is no influence from signal (I chose 0-100 GeV)
                #collecting the new data
                '''
            
                yval, edge = np.histogram(values, range=[100,160], bins=35) #keeping the same bin density
                centre = edge[1:] - 6/7 
            
            #create xvalues to plot the test graphs
                xval = np.linspace(10,100,100)
            
            #using chi squared method to find best fit
            #outer loop is to work on the lambda values (parse through between 25 and 35 as after some trialing we know the value 
            # is between them)
            #inner loop is for the A value, we will trial between 55000 and 60000 as ater testing the value is between those 
            #for i in range(100,160):
                
                y2=stom.get_SB_expectation(bin_centre, 56560, 30, 125, 1.5, 700)
               
                chis = []
                min_chi = []
                min_chi_position = []
                sum = 0
                j = 28
                #y3=(700/(np.sqrt(2.*np.pi)*1.5)*np.exp(-np.power((centre[k] - 125)/1.5, 2.)/2))
                #print('y3',y3)
                while j <= 32: 
                    for i in range(45000, 70000, 10):
                        for k in range(0, len(centre)):
                            psquared = (float((yval[k] - y2[k])/np.sqrt(yval[k]))**2)
                            sum = sum + psquared
                        chis.append(sum)
                        sum = 0
                    
                        #print('chis',chis)
                        #print(psquared)
                    
                    min_chi.append(min(chis))
                    min_chi_position.append(chis.index(min(chis)))
                   
                    #chis = []
                    j = j + 0.1
                    
                  
                    
                smallest_chi_position = min_chi.index(min(min_chi))   
                lamda_value = 28 + smallest_chi_position * 0.1
                A_value = 45000 + min_chi_position[smallest_chi_position] * 10
                
                
                
                
                '''
                #Plotting new test fit data graph
                plt.scatter(centre, yval, color = 'black', s = 25, label = 'Data')
                plt.plot(xval, get_B_expectation(xval, A_value, lamda_value))
                plt.title('Data To find Coefficients')
                yerr = yval*0.01
                '''
                
                #Plots of data and fits
                plt.scatter(bin_centre, bin_heights, color = 'black', s = 25, label = 'Data')
                plt.errorbar(bin_centre, bin_heights, xerr = np.full(shape = 35, fill_value=6/7), yerr = np.sqrt(bin_heights) , ls = 'None', elinewidth = 1, ecolor = 'black', capsize=2)
                plt.plot(bin_centre, y2, label = 'Background Fit')
                
                '''
                #Q2: How does our graph now compare with figure 1?
                Our background fits are very similar, the data tends away from the background fit around the 125 mark, and joins it back again 
                after. There are similar anomolies between the graphs arounf the 140 mark where the number of photons on the bins has a darastic
                decrease
                '''
                
                
                #Formatting of Graph
                plt.xlabel('Invariant Mass/GeV')
                plt.ylabel('Number of Entries')
                plt.legend()
                #plt.legend([l], loc=1)
                plt.show()
                
                '''
                #Q3- I check the goodness of fit by taking the chi squared as a random variable with a PDF that is dependent only on the number of degrees of freedom in the fit.
                minimumchi = min(min_chi) #take the minimum chi value from the min_chi array that was appended in the loops above
                print(minimumchi)
                '''
                
                #Here I calculate the degrees of freedom
                NDOF = len(centre)-2 #Number of degrees of freedom is defined as: number of data points-number of parameters being estimated in the fit. There are 2 fitted parameters here: A and lambda
                print(NDOF)
                
                
                from scipy.stats import chi2 #The scipy chi2 package is used to calculate and plot the PDF of the chi value
                mean, var, skew, kurt = chi2.stats(NDOF, moments='mvsk') #The first 4 moments of the PDF are calculated for the given number of degrees of freedom
                x = np.linspace(chi2.ppf(0.0000000001, NDOF),
                                chi2.ppf(0.9999999999, NDOF), 100) #The PDF has area 1 and is mapped by 100 data points
                plt.plot(x, chi2.pdf(x, NDOF), #Plotting the PDF for the chi value
                       'r-', lw=1, alpha=0.6, label='chi2 pdf') #https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.chi2.html
                plt.ylabel('P(x^2)')
                plt.xlabel('X^2')
                #plt.legend([l], loc=1)
                plt.show()
                from scipy import integrate #This package is used to integrate the PDF between the minimum chi value of the data and 0
                y = chi2.pdf(x, NDOF)
                y_int = integrate.cumulative_trapezoid(y, x, initial=0) #
                def find_nearest(array, value):   #Define a function which finds the closest value to the minimum chi in the array 'x' of 100 points mapping the PDF x-axis
                    array = np.asarray(array)
                    idx = (np.abs(array - value)).argmin()
                    return array[idx]   
                valueinx = find_nearest(x, min(min_chi))
                indexofvalueinx = np.where(x == valueinx)[0]  #Takes the index in the array 'x' of that closest value to the minimum chi
                chiprobability = y_int[indexofvalueinx[0]] #Uses the closest value in 'x' to the minimum chi as bounds for the integration
                print(chiprobability) #Prints the probability that that minimum chi is less than or equal to the calculated minimum chi for the given no. degrees of freedom
                if chiprobability <= 0.9: #Carries out a single-tailed hypothesis test t0 the 10% significance level... higher minimmum chis indicate bad goodness of fit so we only need to test on the upper end
                    print('The goodness of fit is good; we accept the hypothesis that the data agrees with the function chosen for the fit')
                else:
                    print('The goodness of fit is bad; we reject the hypothesis that the data agrees with the function chosen for the fit')
                #if chiprobability==0.05:
                 #  print('correct value', [l])
                  # print('correct value', [l])
                   #print('correct value', [l])
                   #print('corr  ect value', [l])
                   #print('correct value', [l])
                #else:
                 #  continue
                 
                
                # total_chi=list()
                total_chi.append(min(min_chi))
                total_p.append(chiprobability)
                chi3=[108.75543524028593, 121.80180103571264, 125.82370337153803, 123.08144954101137, 116.28469444422147, 108.96726282464401, 107.45602735674696, 112.20334126629784, 113.42409270167765, 104.55211778042404, 92.57499231408603, 88.37810086286204, 96.32866680209942, 111.3772975160466, 123.0842347606455, 124.2657681813777, 116.281139501993, 107.42490311566803, 104.45162615258087, 107.31325917189143, 110.6759763626169, 106.96417572338098, 91.68616354195083, 67.37664044152639, 43.499728676602594, 33.58102542983398, 45.9260046379782, 77.6028961347923, 115.40437465200775, 143.51752769221784, 157.28716289099648, 160.68579625377586, 161.9198171922967, 166.9000994536688, 177.38637386650936, 188.07766332249454, 190.89533372456566, 187.33650613433025, 183.2197574055761, 183.9974448129437, 188.58520473094575, 192.58591236538186, 197.79692622543018, 206.1322361831117, 214.73826336171464, 217.51115806745534, 211.4294611889296, 196.92428968876283, 180.2542021359546, 179.4210869346325, 203.83895919991784, 240.8743964137004, 269.0360036548628, 271.3238324056486, 250.7832264502591, 221.18352831443966, 201.6327131556335, 202.05261105685406, 218.3295131715471, 221.3743184589677]
                
                e=np.linspace(0,np.size(chi3),np.size(chi3))
                plt.scatter(e,chi3)
                plt.show()
                
        
        
                n=n+1
print(total_chi)
print(total_p)
chi3=[108.75543524028593, 121.80180103571264, 125.82370337153803, 123.08144954101137, 116.28469444422147, 108.96726282464401, 107.45602735674696, 112.20334126629784, 113.42409270167765, 104.55211778042404, 92.57499231408603, 88.37810086286204, 96.32866680209942, 111.3772975160466, 123.0842347606455, 124.2657681813777, 116.281139501993, 107.42490311566803, 104.45162615258087, 107.31325917189143, 110.6759763626169, 106.96417572338098, 91.68616354195083, 67.37664044152639, 43.499728676602594, 33.58102542983398, 45.9260046379782, 77.6028961347923, 115.40437465200775, 143.51752769221784, 157.28716289099648, 160.68579625377586, 161.9198171922967, 166.9000994536688, 177.38637386650936, 188.07766332249454, 190.89533372456566, 187.33650613433025, 183.2197574055761, 183.9974448129437, 188.58520473094575, 192.58591236538186, 197.79692622543018, 206.1322361831117, 214.73826336171464, 217.51115806745534, 211.4294611889296, 196.92428968876283, 180.2542021359546, 179.4210869346325, 203.83895919991784, 240.8743964137004, 269.0360036548628, 271.3238324056486, 250.7832264502591, 221.18352831443966, 201.6327131556335, 202.05261105685406, 218.3295131715471, 221.3743184589677]

e=np.linspace(0,np.size(chi3),np.size(chi3))
plt.scatter(e,chi3)
plt.xlabel('X^2 values')
plt.ylabel('Assumed mean Higgs-Boson Rest Mass/GeV')
plt.show()

p_vals=[0.999999999280671, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.999999999280671, 0.9999999990335975, 0.9999999996343134, 0.9999999996847465, 0.999999997325827, 0.9999998340385224, 0.999999277582874, 0.9999999629866357, 0.9999999995591082, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999990335975, 0.999999997325827, 0.9999999990335975, 0.9999999995591082, 0.9999999986674692, 0.9999997595430404, 0.999554382984537, 0.8978500017865763, 0.5662119602551451, 0.9320877427289525, 0.9999827528593399, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465, 0.9999999996847465]
j=np.linspace(100,160,60)

plt.scatter(j,p_vals)
plt.xlabel('Assumed Rest Mass Energy/Gev')
plt.ylabel('P-value')
plt.show()
#j=np.linspace(1,np.size(total_chi),np.size(total_chi))
#plt.scatter(j,total_chi)
#plt.ylabel('X^2')
#plt.xlabel('Trial Number')


