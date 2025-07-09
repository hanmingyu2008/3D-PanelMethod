from math import pi,cos,sin
a = 0.47465
b = 0.31761
c = 0.00124
d = -0.11465
e = -0.4725
f = -0.47668
g = -0.48546
h = 0.0278
h1 = 0.0548
h2 = 0.104
m = 4
theta = 2*pi/(8*m)
x_list = []
r_list = []
for i in range(1,10):
    x = (a*(10-i)+b*i)/10
    x_list.append(x)
    r_list.append(h*i/10)
for i in range(0,11):
    x = (b*(10-i)+c*i)/10
    x_list.append(x)
    r_list.append(h)
for i in range(0,3):
    x = (c*(2-i)+d*i)/2
    x_list.append(x)
    r_list.append(h1)
for i in range(0,11):
    x = (d*(10-i)+e*i)/10
    x_list.append(x)
    r_list.append(h)
for i in range(0,2):
    x = (e*(1-i)+f*i)
    x_list.append(x)
    r_list.append(h2)
for i in range(0,6):
    x = (f*(5-i)+g*i)/5
    x_list.append(x)
    r_list.append(h)
