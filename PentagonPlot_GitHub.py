#plotting multispectral samples from Sentinel-2 images in pentagon plot

#The input is a text file created by extraction of reflectance values fron a Sentinal-2 image in ESA SNAP software



#classes can be shown by markers and colors (input as collumn in the input table)

import math

import matplotlib.pyplot as plt

from collections import OrderedDict

##--------------------------------------------------------------------------input section START

inFile = r'D:\S2_pentagon_data\LagodiGarda\pixEx_S2_MSI_Level-2A_measurements_ALL.txt'


bandPos = [12, 16, 10, 11, 20]  #position of band reflectances for the standard textfile output from PixEx (SNAP)

bandNames = ['R(B4)', 'NIR\n(B8)', 'B(B2)', 'G(B3)', 'SWIR\n(B11)']

posNum = 1     #collumn position with point numbers (integer)



classified = 1        #1: classified, 0: not-classified

posClas = 0           #position of the classification column in the input data

classNames = ['Forest', 'Non-forest vegetation', 'Bare fields', 'Rocks', 'Lake', 'River', 'Snow', 'Clouds', 'Roofs-terracotta', 'Roofs-modern']

colors = ['darkgreen', 'limegreen', 'gold', 'slategrey', 'mediumblue', 'turquoise',  'deepskyblue', 'blueviolet', 'crimson', 'violet']



annotPoints = 0     #annotation of points by pointNum


doAtanStretch = 1

st = 10    #steepness of the curve, in range (0, 5)

sh = -0.4   #horizontal shift of the curve, in range (-0.5, 0.5), negatve value: left, positive: right


gamma = 0.3

if doAtanStretch == 1: gamma = 1    


##--------------------------------------------------------------------------input section END

pi = math.pi

c = [0, 2*pi/5, 4*pi/5, 6*pi/5, 8*pi/5, 0]  #angles of the corner points in pentagon


def atanStretch(x, sh, st):

    #stretch of one value using atan stretch

    y = (math.atan(st *(x - 0.5 - sh)) + (pi/2)) / pi

    offset = (math.atan(st *(0 - 0.5 - sh)) + (pi/2)) / pi

    gain = ((math.atan(st *(1 - 0.5 - sh)) + (pi/2)) / pi) - offset

    y2 = y - offset

    y3 = y2/gain

    return y3


#PLOT HEXAGON

#Corner points of the hexagon with [0,0] in its center and distance corner-center equal 1

xHex = [math.sin(c[0]), math.sin(c[1]), math.sin(c[2]), math.sin(c[3]), math.sin(c[4]), math.sin(c[0])]
yHex = [math.cos(c[0]), math.cos(c[1]), math.cos(c[2]), math.cos(c[3]), math.cos(c[4]), math.cos(c[0])]

plt.plot(xHex, yHex, 'k-', linewidth = 0.5)

#Lables for the pentagon corner points

fs = 10
ad = 1.02

plt.annotate(bandNames[0], (xHex[0]*ad, yHex[0]*ad), fontsize=fs, ha='center')
plt.annotate(bandNames[1], (xHex[1]*ad, yHex[1]*ad), fontsize=fs, ha='left')
plt.annotate(bandNames[2], (xHex[2]*ad, yHex[2]*ad), fontsize=fs, va='top', ha='center')
plt.annotate(bandNames[3], (xHex[3]*ad, yHex[3]*ad), fontsize=fs, va='top', ha='center')
plt.annotate(bandNames[4], (xHex[4]*ad, yHex[4]*ad), fontsize=fs, ha='right')



#plot star lines

for i in range(0,6):

    plt.plot([0,xHex[i]], [0 ,yHex[i]], 'k', linewidth = 0.1)


#hexagon lines for selected reflectance values

#rList = [0.25, 0.5, 0.75, 1.0]
#anLst = ['0.25', '0.5', '0.75', '1.0']

rList = [0.1, 0.2, 0.4, 1.0]
anLst = ['0.1', '0.2', '0.4', '1.0']

annotX = []
annotY = []
ii = 0

xH = [0,0,0,0,0,0]
yH = [0,0,0,0,0,0]

for r in rList:

    for i in range(0,6):

        #inner pentagon lines

        if doAtanStretch == 1:

            r2 = atanStretch(r, sh, st)

        else:

            r2 = r**gamma

        xH[i] = math.sin(c[i]) * r2
        yH[i] = math.cos(c[i]) * r2


    plt.plot(xH, yH, 'k-', linewidth = 0.1)


    #annotations of the inner pentagon lines (0.25, 0.5, 0.75, 1.0)

    adj = 0.02

    angList = [pi/5, 3*pi/5, 5*pi/5, 7*pi/5, 9*pi/5]

    rotLab = [-40, 72, 0, -72, 36]

    sel = 3          #position of annotations: 0: upper right, 1:bottom right, 2: bottom, 3: bottom left


    ra = r * math.cos(pi/5)

    fi = angList[sel]

    p5 = 2*pi/5

    ff = (fi % p5)   #modulo

    r0 = ra * math.cos(pi/5 - ff) / math.cos(pi/5)


    if doAtanStretch == 1:

        r2 = atanStretch(r0, sh, st)

        v = r2 * ra / r0 + adj

    else:

        v = (r0**gamma) * ra / r0 + adj



    annotX.append(v*math.sin(angList[sel]))

    annotY.append(v*math.cos(angList[sel]))


    plt.annotate(anLst[ii], (annotX[ii], annotY[ii]), rotation=rotLab[sel], color = 'black', size=7, va='center')

    ii = ii+1



#PROCESS SAMPLE DATA

f = open(inFile, 'r')

ss = [0,0,0,0,0]

i = 0

for lin in f:

    if i > 6:

        #read spectral values and point number form a text file

        for j in range(0,5):

            pointNum = int(lin.split("\t")[posNum])

            if lin.split("\t")[posClas]:
                clas = int(lin.split("\t")[posClas])
                
            else:
                clas = -1

            bp = bandPos[j]

            if lin.split("\t")[bp]:
                ss[j] = float(lin.split("\t")[bp])
                
            else:
                ss[j] = 999
                print('input value missing, replaced by 999. Point number, band :', pointNum, bandNames[j] )



        if ss[0] != 999 and ss[1] != 999 and ss[2] != 999 and ss[3] != 999 and ss[4] != 999:


            #normalization

            s = [0,0,0,0,0]

            for ii in range(0,5):

                 s[ii] = ss[ii] / (ss[0] + ss[1] + ss[2] + ss[3] + ss[4])


            #vector addition

            xTot = s[0] * math.sin(c[0]) + s[1] * math.sin(c[1]) + s[2] * math.sin(c[2]) + s[3] * math.sin(c[3]) + s[4] * math.sin(c[4])

            yTot = s[0] * math.cos(c[0]) + s[1] * math.cos(c[1]) + s[2] * math.cos(c[2]) + s[3] * math.cos(c[3]) + s[4] * math.cos(c[4])


            #gamma scaling or arkustangens stretch of the resulting radius (conversion to polar coordinates first)

            r = (xTot**2  + yTot**2)**0.5

            fi = math.atan2(xTot, yTot)

            p5 = 2*pi/5

            ff = (fi % p5)   #modulo

            r0 = r * math.cos(pi/5 - ff) / math.cos(pi/5)

            if doAtanStretch == 1:

                r1 = atanStretch(r0, sh, st)    #tangens stretch

                r2 = r1 * r / r0

            else:

                r2 = (r0**gamma) * r / r0   #gamma scaling


            xTot = r2 * math.sin(fi)

            yTot = r2 * math.cos(fi)


            #plot points

            if classified == 1:

                plt.plot(xTot, yTot,  '.', color = colors[clas-1], label=classNames[clas-1])

            else:

                plt.plot(xTot, yTot, 'k+')
                
                

            if annotPoints == 1: plt.annotate(pointNum, (xTot, yTot+0.04))


    i = i+1
    

#finish plotting

plt.gca().set_aspect('equal')

plt.axis('off')


handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
leg = plt.legend(by_label.values(), by_label.keys(), loc = 'lower left', bbox_to_anchor=(1.05, 0.0))

leg.set_draggable(state=True)


plt.show()

print('*** finished ***')
