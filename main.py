import scipy.optimize as optimize
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
from pylab import loadtxt
import math
import statistics as stat

def get_counts(filename):
      # This function unpacts the data from the text files and computs
  # rstep: the distance between the last point and the current point
  # dis: the total sqaured displacement of the particle


  data=loadtxt(filename, usecols=(0,1), skiprows=2, unpack=True)
  xdata=data[0]
  ydata=data[1]

  # Initializing reference origin

  
  rstep = []
  dis = [0]

  # Converting from pixles to micrometers
  for i in range (len(xdata)): 
    xdata[i] = xdata[i] * 0.1155
    #0.12048
    ydata[i] = ydata[i] * 0.1155

  x = xdata[0]
  y = ydata[0]
  

  for i in range (1,len(xdata)-1):
      #Initializing changing origin
      rstep.append(np.sqrt((xdata[i] - xdata[i-1])**2 + (ydata[i] - ydata[i-1])**2))
      dis.append((xdata[i] - x)**2 +(ydata[i] - y)**2)
  return rstep, dis

def get_data():

    step2, dis2 = get_counts('run2_data.txt') 
    step3, dis3 = get_counts('run3_data.txt') 
    step4, dis4 = get_counts('run4_data.txt') 
    step5, dis5 = get_counts('run5_data.txt') 
    step6, dis6 = get_counts('run6_data.txt') 
    step7, dis7 = get_counts('run7_data.txt') 
    step9, dis9 = get_counts('run9_data.txt') 
    step11, dis11 = get_counts('run11_data.txt') 
    step12, dis12 = get_counts('run12_data.txt') 
    step13, dis13 = get_counts('run13_data.txt') 

    all_step = []
    ave_step = []
    deviation = []

    for i in range(0,119):
        sum = dis2[i] + dis3[i] + dis4[i] + dis5[i] + dis6[i] + dis7[i] + dis9[i] + dis11[i] + dis12[i] +dis13[i]
        ave_step.append(sum/10)
        deviation.append(stat.stdev([dis2[i], dis3[i], dis4[i], dis5[i], dis6[i], dis7[i], dis9[i], dis11[i], dis12[i], dis13[i]]))

    for i in range(0, 118):
        all_step.append(step2[i])
        all_step.append(step3[i])
        all_step.append(step4[i])
        all_step.append(step5[i])
        all_step.append(step6[i])
        all_step.append(step7[i])
        all_step.append(step9[i])
        all_step.append(step11[i])
        all_step.append(step12[i])
        all_step.append(step13[i])
    
    return all_step, ave_step, deviation

def my_func(x, a):
  # Function for the linear fit
  return 4*a*x

def my_func2(x, D):
  # Function for the probability density fit
  return (np.exp((-x**2)/(4*D*0.5))*x/(2*D*0.5))

def lin_fit(ave_step):
  # Linear Fit function
  # Returns time (x_values), residuals, and y_values
  # Also calcualtes the slope and its error
  time = [0]
  new_y =[]
  for i in range (1, 119):
    time.append(0.5*i)

  popt, pcov = optimize.curve_fit(my_func, time, ave_step)
  m = popt[0]
  m_err = np.sqrt(pcov[0][0])
  print("Linear Fit D: ", m)
  print("Linear Fit D Error: ", m_err)
  res = []
  for i in range (0, 119):
    new_y.append(4*m*time[i])
    new_y[i] = new_y[i]
    res.append(new_y[i] - ave_step[i])
  return time, res, new_y

def hist_fit (bins, n):
  # Histogram fit function
  popt, pcov = optimize.curve_fit(my_func2, bins, n)
  D = popt[0]
  D_err = np.sqrt(pcov[0][0])
  return D, D_err

def Est2dt (all_step):
  N = len(all_step)
  sum = 0
  for i in range(0, N):
    sum = sum + all_step[i]**2
  
  sum = sum/N
  sum = sum/2
  return sum

def chi_msd(ave_step, new_y, yerr):
  # need to store amount for final sum
  elements = []
  chi_1 = 0
  sum = 0
  for i in range(0, len(ave_step)):
    numerator = (ave_step[i] - new_y[i])**2
    denominator = yerr[i]**2
    sum = sum + numerator;
    elements.append(numerator/denominator)
  for n in range(0, len(elements)):
    chi_1 = chi_1 + elements[n]
  return chi_1

def plotting(time, res, new_y, all_step, distibution):
  # Plotting function
  yerr = []
  xerr = []
  for i in range(0, len(res)):
    yerr.append(0.1)
    xerr.append(0.03)
  binnum = 100
  fig, (p1, p2) = plt.subplots(2, 1)
  fig.subplots_adjust(hspace=0.6)
  fig.set_figheight(30)
  fig.set_figwidth(15)

  fig2, (p3, p4) = plt.subplots(2, 1)
  fig2.subplots_adjust(hspace=0.6)
  fig2.set_figheight(30)
  fig2.set_figwidth(15)
  n, bins, patches = p3.hist(all_step, density = True , bins = binnum)
  histerr = []
  for i in range(0, len(n)):
    histerr.append((n[i]*len(all_step))**(1/2)/len(all_step))
    if histerr[i] == 0:
        histerr[i] = 0.001
    
  mid = 0.5*(bins[1:] + bins[:-1])
  p3.errorbar(mid, n, histerr, fmt='none', ecolor = "blue")
  #p1.errorbar(time, ave_step,yerr, xerr, fmt="o")
  p1.errorbar(time, ave_step,distibution, xerr, fmt="o")
  p1.plot(time, new_y)
  p1.title.set_text("Linear Fit Mean Squared Distance vs Time")
  #p2.errorbar(time, res, yerr, xerr, fmt = 'o')
  p2.errorbar(time, res, distibution, xerr, fmt = 'o')
  p2.title.set_text("Residuals of Linear Fit")
  p3.set_ylabel("Probaility of Occurances")
  p3.title.set_text("Probability Density of Step Sizes")
  p3.set_xlabel("Step Size [Micro Meters]")
  p1.set_xlabel("Time [s]")
  p1.set_ylabel("Mean Squared Distance [Micro Meters]")
  p1.legend("Scatter \n Fitted curve = {a}x + {b}")
  p2.set_ylabel("Residuals of Linear Fit [Micro Meters]")
  p2.set_xlabel("Time[s]")
  p4.title.set_text("Residuals of Gaussian Fit")
  
  all_step, ave, distribution = get_data()
  sum = Est2dt(all_step)
  print("Maximum Likihood D: ", sum)
  
  # Changing bins to average step of each bin
  x_val = []
  for i in range(0,binnum):
   x_val.append((bins[i]+ bins[i + 1])/2)
  # Points of the histogram


  D, D_err = hist_fit(x_val, n)
  print("Hitogram Fit D: ", D)
  print("Histogram Fit D Error: ", D_err)
  
  # Calculating points of fit
  y_val = []
  x = np.arange(0, x_val[-1], 0.01)
  for i in range(0, len(x)):
    y_val.append(np.exp((-x[i]**2)/(4*sum*0.5))*x[i]/(2*sum*0.5))
  #p3.plot(x, y_val)
  y = my_func2(mid, D);
  p4.errorbar(mid, n-y, histerr, fmt = "o")

  zero = []
  lin = np.arange(0, 70, 10)
  z = [0, 0, 0, 0, 0, 0, 0]
  for i in range(0, binnum):
    zero.append(0)
  p2.plot(lin, z)
  p4.plot(mid, zero)
  p4.set_ylabel("Residuals")
  p4.set_xlabel("Step Size [Micro Meters]")

  plt.show()

  #function to calculate chi of probability
  print("Chi Square Histrogram: ",chi_msd(n, y_val, histerr)/(len(n)-1))


all_step, ave_step, distribution = get_data()
time, res, new_y = lin_fit(ave_step)
plotting(time, res, new_y, all_step, distribution)
sum = Est2dt(all_step)
print("Maximum Likihood D: ", sum)

#different bins fall within eror

#function to calculate chi of mean squared distance graph
yerr = []
for i in range(0, len(new_y)):
  yerr.append(0.1)
  if distribution[i] < 0.1:
      distribution[i] = 0.1
print("Chi Square Linear: ", chi_msd(ave_step, new_y, distribution)/119)




all_step, ave_step, distribution = get_data()
time, res, new_y = lin_fit(ave_step)
plotting(time, res, new_y, all_step, distribution)
sum = Est2dt(all_step)
print("Maximum Likihood D: ", sum)

#different bins fall within eror

#function to calculate chi of mean squared distance graph
yerr = []
for i in range(0, len(new_y)):
  yerr.append(0.1)
  if distribution[i] < 0.1:
      distribution[i] = 0.1;
print("Chi Square Linear: ", chi_msd(ave_step, new_y, distribution)/119)

