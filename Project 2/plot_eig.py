import matplotlib.pyplot as plt
import sys
import numpy as np

def analytical_eigenvectors(innerpoints):
    eigvec = np.zeros((innerpoints,innerpoints))
    for i in range(1,innerpoints+1):
        for j in range(1,innerpoints+1):
            eigvec[i-1,j-1] = np.sin((i*j*np.pi)/(innerpoints+1))
    return eigvec


if len(sys.argv) < 1:
    print("Please write the name of the file you wish to plot")
else:
    filename = str(sys.argv[1])
    myfile = open(filename, "r")
    eigenvectors = []
    for i in myfile:
        i.strip()
        if i != "\n":
            eigenvectors.append(float(i))
    myfile.close()
    n = int(sys.argv[1])
    n_inner = n-1
    eigenvectors1 = eigenvectors[:n_inner]
    eigenvectors2 = eigenvectors[n_inner:-n_inner]
    eigenvectors3 = eigenvectors[-n_inner:]

    # Add end-points
    eigenvectors1.append(0)
    eigenvectors2.append(0)
    eigenvectors3.append(0)
    eigenvectors1.insert(0,0)
    eigenvectors2.insert(0,0)
    eigenvectors3.insert(0,0)

    x = []
    h = 1/n
    for j in range(n+1):
        x.append(j*h)

    eigvec_analytical = analytical_eigenvectors(n_inner)
    norm = np.linalg.norm(eigvec_analytical[0])

    # Normalise analytical eigenvectors
    eigvec1 = eigvec_analytical[0]/norm
    eigvec2 = eigvec_analytical[1]/norm
    eigvec3 = *eigvec_analytical[2]/norm

    plt.plot(x, eigenvectors1, label="1")
    plt.plot(x, eigenvectors2, label="2")
    plt.plot(x, eigenvectors3, label="3")
    plt.plot(x, np.concatenate([[0],eigvec1,[0]]), "--r", label="analytical1")
    plt.plot(x, np.concatenate([[0],eigvec2,[0]]), "--c",label="analytical2")
    plt.plot(x, np.concatenate([[0],eigvec3,[0]]), "--k",label="analytical3")

    plt.xlabel("$x_i$")
    plt.ylabel("$v_i$")
    plt.legend()
    plt.savefig(filename + "minus.pdf")
