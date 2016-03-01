#!/usr/bin/python
__author__ = "morganlnance"



def calc_distance( vec1, vec2 ):
    """
    Calculates the distance between two lists of 3D coordinates.
    :param vec1: list( xyz coordinates defining point 1 )
    :param vec1: list( xyz coordinates defining point 2 )
    :return dist: float( the distance between the two points in 3D space )
    """
    from math import sqrt, pow
    
    # extract the xyz coordinates from the two points
    x1 = vec1[0]
    y1 = vec1[1]
    z1 = vec1[2]
    
    x2 = vec2[0]
    y2 = vec2[1]
    z2 = vec2[2]
    
    # calculate the distance between the two points
    dist = sqrt( pow( x2 - x1, 2 ) + pow( y2 - y1, 2 ) + pow( z2 - z1, 2 ) )
    
    return dist
