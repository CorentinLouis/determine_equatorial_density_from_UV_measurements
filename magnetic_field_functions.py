import numpy
from coordinates_system_transformation import spherical_to_cartesian, cartesian_to_spherical



def read_b_model(file):
    data_file = numpy.loadtxt(file,dtype="str", delimiter= ',')
    r = data_file[:,0].astype(float)
    theta = data_file[:,1].astype(float)
    phi = data_file[:,2].astype(float)
    B_total = data_file[:,3].astype(float)

    return (r, theta, phi, B_total)



def length_magnetic_field_line(coord1, coord2, coord3, spherical = True, cartesian = False):
    """
    If 'spherical = True':
        (coord1, coord2, coord3) = (r, theta, phi) the spherical coordinates
    If 'cartesian = True':
        (coord1, coord2, coord3) = (x,y,z) the cartesian coordinates)
    """
    
    if (spherical == True) and (cartesian == True):
        print('ERROR')

    if spherical:
        (x,y,z) = spherical_to_cartesian(coord1,coord2,coord3)
    else:
        (x,y,z) = (coord1, coord2, coord3)
    
    
    cartesian_points = numpy.zeros((len(coord1), 3))
    for i in range(cartesian_points.shape[0]):
        cartesian_points[i,:] = [x[i], y[i], z[i]]
    length_mfl_B = 0.0

    for i in range(1, len(cartesian_points)):
        length_mfl_B += numpy.linalg.norm(cartesian_points[i] - cartesian_points[i - 1])

    return length_mfl_B




def delta_length_magnetic_field_line(coord1, coord2, coord3, spherical = True, cartesian = False):
    """
    If 'spherical = True':
        (coord1, coord2, coord3) = (r, theta, phi) the spherical coordinates
    If 'cartesian = True':
        (coord1, coord2, coord3) = (x,y,z) the cartesian coordinates)
    """
    
    if (spherical == True) and (cartesian == True):
        print("ERROR")

    if spherical:
        (x,y,z) = spherical_to_cartesian(coord1,coord2,coord3)
    else:
        (x,y,z) = (coord1, coord2, coord3)
    
    
    cartesian_points = numpy.zeros((len(coord1), 3))
    for i in range(cartesian_points.shape[0]):
        cartesian_points[i,:] = [x[i], y[i], z[i]]

    delta_length_mfl_B = numpy.zeros((len(coord1)-1))
    for i in range(0, len(delta_length_mfl_B)):
        delta_length_mfl_B[i] = numpy.linalg.norm(cartesian_points[i+1] - cartesian_points[i])

    return delta_length_mfl_B
