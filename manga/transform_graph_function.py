import math
def transform(sx,sy,thet):
    theta = float(thet)
    x = float(sx)
    y = float(sy)
    xprime = (x * math.cos(math.radians(theta))) -(y * (math.sin(math.radians(theta))))
    yprime = (y * (math.cos(math.radians(theta)))) + (x * (math.sin(math.radians(theta))))
    return(str(round(xprime,2)),str(round(yprime,2)))


data = input('please enter your x-coordinate,y-coordinate,theta: ')
xd, yd, thetad = data.split(',')
print(transform(xd,yd,thetad))
