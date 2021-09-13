import matplotlib.pyplot as plt
import sys
import numpy as np
import glob


filenames = sorted(glob.glob('N?'))
max_rerror = {}

for filename in filenames: # for N1, N2, N3
    #print(filename)
    myfile = open(filename, "r")
    e = []
    x = []
    n = []
    aerror = []
    rerror = []
    for i in myfile:
        i = i.split()
        x.append(i[0])
        n.append(i[1])
        e.append(i[2])
        aerror.append(i[3])
        rerror.append(i[4])
    max_rerror[filename] = max(rerror)
    # Plotting numerical vs exact solutions
    #plt.plot(x,n, label=filename)
    #if filename == filenames[-1]: # Plot exact solution for maximum number of N
        #plt.plot(x,e, label="eksakt")
    #plt.legend()
    #plt.xlabel("x")
    #plt.ylabel("u(x)")
    #plt.savefig(filename + ".pdf")

    # Plotting absolute error
    #plt.plot(x, aerror, label=filename)
    #plt.xlabel("x")
    #plt.ylabel("$log_{10}(|u_i-v_i|)$")
    #plt.legend()
    #plt.savefig("abs_error for" + filename + ".pdf")

    # Plotting relative error
    plt.plot(x, rerror, label=filename)
    plt.xlabel("x")
    plt.ylabel("$log_{10}\epsilon$")
    plt.legend()
    plt.savefig("rel_error for" + filename + ".pdf")

#print(max_rerror)




"""
if len(sys.argv) < 1:
    print("Please write the name of the file you wish to plot")
else:
    filename = str(sys.argv[1])
    myfile = open(filename, "r")
    e = []
    x = []
    n = []
    aerror = []
    rerror = []
    for i in myfile:
        i = i.split()
        x.append(i[0])
        n.append(i[1])
        e.append(i[2])
        aerror.append(i[3])
        rerror.append(i[4])

    # Plotting numerical vs exact solutions
    #plt.plot(x,n, label="numerisk")
    #plt.plot(x,e, label="eksakt")
    #plt.xlabel("x")
    #plt.ylabel("u(x)")
    #plt.savefig("exact.pdf")
    #plt.savefig(str(sys.argv[1]) + ".pdf")
    #plt.legend()

    # Plotting errors
    #plt.plot(x, aerror)
    #plt.xlabel("x")
    #plt.ylabel("$log_{10}(|u_i-v_i|)$")
    #plt.savefig("abs_error" + str(sys.argv[1]) + ".pdf")

    #plt.plot(x,rerror)
    #plt.show()
    myfile.close()
"""
