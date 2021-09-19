import numpy as np
import matplotlib.pyplot as plt


#This first section of code simulates disease spreading using the SIS-model
#Define initial values
#Population size
N = 90000
#Recovery rate of infected individuals
a = 0.125
#Infection rate of susceptible individuals
B = 17*a
#Initial no. Susceptible individuals witin N
I = 1000
#Initial no. Infected individuals witin N
S = N-I
#x is the number of days the program simulates, z is the number of time steps in time period x.
x=30
z=10000
#Create an array of time steps that will be used for plotting later
t = np.linspace(0,x,z)
#Define step size
h = x/z
#This value will be used for the proof that SIS tends to 1-alpha/beta
d=1-(a/B)

#Define lists for h, S and I
Slist=[]
Ilist=[]
flist=[]

#Define the functions for I and S that have been solved with Euler method
def sus(S,I):
    return S-B*S*I*h/N + a*h*I
    
def inf(S,I):
    return I+B*S*I*h/N - a*I*h

#For loop that run the functions the desired amount of times
for i in range(0,z):
     #These consants prevent the program from using value n+1 instead of n each loop
     X=S
     Y=I
     #Run functions and append the result in lists
     S=sus(X, Y)
     Slist.append(S)
     I=inf(X, Y)
     Ilist.append(I)
     #Append the fraction of I.
     flist.append(I/N)

#Plot the lists over time t
plt.figure()
plt.plot(t, Slist, 'r', label='Susceptible')
plt.plot(t, Ilist, 'b', label='Infected')
plt.xlabel('Days')
plt.ylabel('Individuals')
plt.legend(bbox_to_anchor=(1.1, 1.05))
plt.show()

#Check if the fraction of infected goes to 1-α/β by plotting it
plt.figure()
plt.plot(t, flist, 'r', label='Fraction f')
plt.axhline(1-(a/B), label='1-α/β')
plt.ylabel('Time')
plt.legend(bbox_to_anchor=(1.1, 1.05))
plt.show()


#SIR with travel parameters. Travel between Malmö (m) and Lund (n)
#Define the two population sizes.
Nn=90000
Nm=340000
#Define intial values for S, I and R in each population
In=1000
Im=1000
Rn=0
Rm=0
Sn=Nn-In-Rn
Sm=Nm-Im-Rm
#Define vaccination rate constant
y=0.01
#Define rate constants for travel from population n to m
wnm=0.1
#Define rate constants for travel from population m to n
wmn=0.2
#Define new lists for S, I and R
Snlist=[]
Smlist=[]
Inlist=[]
Imlist=[]
Rnlist=[]
Rmlist=[]
#Define constants for error estimation as well
NnE=Nn
NmE=Nm
InE=In
ImE=Im
RnE=Rn
RmE=Rm
SnE=Sn
SmE=Sm

#Define functions for SIR with travel.
def susT(N1,S1,S2,I1,wmn,wnm):
    return S1-B*S1*I1*h/N1 - y*h*S1 + h*(wmn*S2-wnm*S1)
def infT(N1,S1,I1,I2,wmn,wnm):
    return I1+B*S1*I1*h/N1 - a*I1*h + h*(wmn*I2-wnm*I1)
def recT(S1,I1,R1,R2,wmn,wnm):
    return R1 + h*y*S1 + h*a*I1 + h*(wmn*R2-wnm*R1)

#For Loop that runs the functions z times
for i in range(0,z):
       # Set the populations to be the sum of S, I and R.
       Nn=Sn+In+Rn
       Nm=Sm+Im+Rm
       # These constants hold the values for S, I and R so that the functions don't use n+1 values until the whole loop is done.
       Xn=Sn
       Xm=Sm
       Yn=In
       Ym=Im
       Zn=Rn
       Zm=Rm
       # Run the recursion functions for S, I and R in each population and append the values.
       Sn=susT(Nn,Xn,Xm,Yn,wmn,wnm)
       Sm=susT(Nm,Xm,Xn,Ym,wnm,wmn)
       Snlist.append(Sn)
       Smlist.append(Sm)
       In=infT(Nn,Xn,Yn,Ym,wmn,wnm)
       Im=infT(Nm,Xm,Ym,Yn,wnm,wmn)
       Inlist.append(In)
       Imlist.append(Im)
       Rn=recT(Xn,Yn,Zn,Zm,wmn,wnm)
       Rm=recT(Xm,Ym,Zm,Zn,wnm,wmn)
       Rnlist.append(Rn)
       Rmlist.append(Rm)
      
#Plot the lists over t.
plt.figure()
plt.plot(t, Snlist, 'r', label='Susceptible in Lund')
plt.plot(t, Inlist, 'b', label='Infected in Lund')
plt.plot(t, Rnlist, 'g', label='Recovered in Lund')
plt.plot(t, Smlist, 'r',linestyle='dashed', label='Susceptible in Malmö')
plt.plot(t, Imlist, 'b',linestyle='dashed', label='Infected in Malmö')
plt.plot(t, Rmlist, 'g',linestyle='dashed', label='Recovered in Malmö')
plt.xlabel('Days')
plt.ylabel('Individuals')
plt.legend(bbox_to_anchor=(1.1, 1.05))
plt.show()

#Error estimation

#For error estimation, we run the code a second time, but with h/2
h=h/2

Snlist2=[]
Smlist2=[]
Inlist2=[]
Imlist2=[]
Rnlist2=[]
Rmlist2=[]
#Define lists that hold the error values.
SnlistE=[]
SmlistE=[]
InlistE=[]
ImlistE=[]
RnlistE=[]
RmlistE=[]

for i in range(0,2*z):
      #Set the populations to be the sum of S, I and R
      NnE=SnE+InE+RnE
      NmE=SmE+ImE+RmE
      #These constants hold the values for S, I and R so that the functions don't use n+1 values until the whole loop is done.
      XnE=SnE
      XmE=SmE
      YnE=InE
      YmE=ImE
      ZnE=RnE
      ZmE=RmE
      #Run the recursion functions with h/2 and append the difference betweeen values with h and values with h/2 to get the error
      #If statement makes the program only append every other value
      SnE=susT(NnE,XnE,XmE,YnE,wmn,wnm)
      SmE=susT(NmE,XmE,XnE,YmE,wnm,wmn)
      InE=infT(NnE,XnE,YnE,YmE,wmn,wnm)
      ImE=infT(NmE,XmE,YmE,YnE,wnm,wmn)
      RnE=recT(XnE,YnE,ZnE,ZmE,wmn,wnm)
      RmE=recT(XmE,YmE,ZmE,ZnE,wnm,wmn)
      #We only append every other values so that we can compute the error in each whole h step
      if(i%2)==0:
          pass
      else:
          Snlist2.append(SnE)
          Smlist2.append(SmE)
          Inlist2.append(InE)
          Imlist2.append(ImE)
          Rnlist2.append(RnE)
          Rmlist2.append(RmE)
          
#Subract the values using h with the corresponding values using h/2. We then append these in a new list.
errorSn=zip(Snlist2,Snlist)
for Snlist2_i, Snlist_i in errorSn:
    SnlistE.append(Snlist2_i-Snlist_i)

errorSm=zip(Smlist2,Smlist)
for Smlist2_i, Smlist_i in errorSm:
    SmlistE.append(Smlist2_i-Smlist_i) 

errorIn=zip(Inlist2,Inlist)
for Inlist2_i, Inlist_i in errorIn:
    InlistE.append(Inlist2_i-Inlist_i) 

errorIm=zip(Imlist2,Imlist)
for Imlist2_i, Imlist_i in errorIm:
    ImlistE.append(Imlist2_i-Imlist_i) 

errorRn=zip(Rnlist2,Rnlist)
for Rnlist2_i, Rnlist_i in errorRn:
    RnlistE.append(Rnlist2_i-Rnlist_i) 

errorRm=zip(Rmlist2,Rmlist)
for Rmlist2_i, Rmlist_i in errorRm:
    RmlistE.append(Rmlist2_i-Rmlist_i)


       
#Plot the error lists over time
plt.figure()
plt.plot(t, SnlistE, 'r', label='S error in Lund')
plt.plot(t, InlistE, 'b', label='I error in Lund')
plt.plot(t, RnlistE, 'g', label='R error in Lund')
plt.plot(t, SmlistE, 'r',linestyle='dashed', label='S error in Malmö')
plt.plot(t, ImlistE, 'b',linestyle='dashed', label='I error in Malmö')
plt.plot(t, RmlistE, 'g',linestyle='dashed', label='R error in Malmö')
plt.xlabel('Days')
plt.ylabel('Error')
plt.legend(bbox_to_anchor=(1.1, 1.05))
plt.show()