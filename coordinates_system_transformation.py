from numpy import pi, cos, sin, arccos, arctan2, sqrt

def spherical_to_cartesian(r,theta,phi):
    phi = (360 - phi) /180 * pi
    theta = theta / 180 * pi
    x = r * sin(theta) * cos(phi)
    y = r * sin(theta) * sin(phi)
    z = r* cos(theta)
    return (x,y,z)


def cartesian_to_spherical(x,y,z):
    r = sqrt(x**2+y**2+z**2)

    theta = arccos(z / r)  # theta is between 0 and pi

    phi = arctan2(y, x)  # phi is between -pi and pi
    phi[phi < 0] += 2 * pi
    phi = 360-(phi/pi*180) #(go to SystIII RH longitude)    
    theta = theta/pi*180
    return r, theta, phi