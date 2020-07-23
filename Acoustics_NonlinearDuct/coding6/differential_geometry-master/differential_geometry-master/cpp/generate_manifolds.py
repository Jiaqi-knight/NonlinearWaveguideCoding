#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, getopt, re, os, time, csv, codecs
from datetime import datetime
from random import randint
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def write_csv(dim,manifold,outdir,outfile):
    
    if outdir[-1] == '/' and outfile[0] != '/' and '.' not in outfile:
        outfilepath = outdir + outfile + '.csv'
    elif outdir[-1] != '/' and outfile[0] == '/' and '.' not in outfile:
        outfilepath = outdir + outfile + '.csv'
    elif outdir[-1] != '/' and outfile[0] != '/' and '.' not in outfile:
        outfilepath = outdir + '/' + outfile + '.csv'
    elif outdir[-1] == '/' and outfile[0] == '/' and '.' not in outfile:
        outfilepath = outdir + outfile[1:] + '.csv'
    elif outdir[-1] == '/' and outfile[0] != '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + newoutfile + '.csv'
    elif outdir[-1] != '/' and outfile[0] == '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + newoutfile + '.csv'
    elif outdir[-1] != '/' and outfile[0] != '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + '/' + newoutfile + '.csv'
    elif outdir[-1] == '/' and outfile[0] == '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + newoutfile[1:] + '.csv'
    
    with open(outfilepath,'w') as CSVfile:
        for element in manifold:
            for index, coordinate in enumerate(element):
                CSVfile.write(str(coordinate))
                if index == dim-1:
                    CSVfile.write('\n')
                else:
                    CSVfile.write(',')

def write_inp(params,outdir,outfile):
    
    if outdir[-1] == '/' and outfile[0] != '/' and '.' not in outfile:
        outfilepath = outdir + outfile + '.inp'
    elif outdir[-1] != '/' and outfile[0] == '/' and '.' not in outfile:
        outfilepath = outdir + outfile + '.inp'
    elif outdir[-1] != '/' and outfile[0] != '/' and '.' not in outfile:
        outfilepath = outdir + '/' + outfile + '.inp'
    elif outdir[-1] == '/' and outfile[0] == '/' and '.' not in outfile:
        outfilepath = outdir + outfile[1:] + '.inp'
    elif outdir[-1] == '/' and outfile[0] != '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + newoutfile + '.inp'
    elif outdir[-1] != '/' and outfile[0] == '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + newoutfile + '.inp'
    elif outdir[-1] != '/' and outfile[0] != '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + '/' + newoutfile + '.inp'
    elif outdir[-1] == '/' and outfile[0] == '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + newoutfile[1:] + '.inp'
    
    with open(outfilepath,'w') as INPfile:
        INPfile.write('// Datafile features\n')
        INPfile.write('//\n')
        INPfile.write('  D = ' + str(params[0]) + '                                       // Space dimensions\n')
        INPfile.write('  Nx = ' + str(params[1]) + '                                    // Number of points in x-direction\n')
        INPfile.write('  Ny = ' + str(params[2]) + '                                    // Number of points in y-direction\n')
        INPfile.write('  Nz = ' + str(params[3]) + '                                    // Number of points in z-direction\n')
        INPfile.write('//\n')
        INPfile.write('// Input filenames\n')
        INPfile.write('//\n')
        INPfile.write('  EndMapManFileName = ' + str(params[4]) + '      // CSV file with coordinates of points describing the co-domain manifold\n')
        INPfile.write('  BaseMapManFileName = ' + str(params[5]) + '    // CSV file with coordinates of points describing the domain manifold\n')
        INPfile.write('//\n')
        INPfile.write('// Output filepath')
        INPfile.write('//')
        INPfile.write('   OutputDirectory = ' + str(params[6]) + '                        // Relative path (to working directory) of output directory')
        INPfile.write('   OutputFileName = ' + str(params[7]) + '              // Filename to save results')
        INPfile.write('//')

def square(l,cx,cy,Nx,Ny):
    coords = np.zeros((Nx*Ny,2))
    for j in range(Ny):
        for i in range(Nx):
            coords[i+j*Nx,0] = cx-0.5*l+i*l/(Nx-1)
            coords[i+j*Nx,1] = cy-0.5*l+j*l/(Ny-1)
    return coords

def rectangle(lx,ly,cx,cy,Nx,Ny):
    coords = np.zeros((Nx*Ny,2))
    for j in range(Ny):
        for i in range(Nx):
            coords[i+j*Nx,0] = cx-0.5*lx+i*lx/(Nx-1)
            coords[i+j*Nx,1] = cy-0.5*ly+j*ly/(Ny-1)
    return coords

def cube(l,cx,cy,cz,Nx,Ny,Nz):
    coords = np.zeros((Nx*Ny*Nz,3))
    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                coords[i+j*Nx+k*Nx*Ny,0] = cx-0.5*l+i*l/(Nx-1)
                coords[i+j*Nx+k*Nx*Ny,1] = cy-0.5*l+j*l/(Ny-1)
                coords[i+j*Nx+k*Nx*Ny,2] = cz-0.5*l+k*l/(Nz-1)
    return coords

def parallelepiped(lx,ly,lz,cx,cy,cz,Nx,Ny,Nz):
    coords = np.zeros((Nx*Ny*Nz,3))
    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                coords[i+j*Nx+k*Nx*Ny,0] = cx-0.5*lx+i*lx/(Nx-1)
                coords[i+j*Nx+k*Nx*Ny,1] = cy-0.5*ly+j*ly/(Ny-1)
                coords[i+j*Nx+k*Nx*Ny,2] = cz-0.5*lz+k*lz/(Nz-1)
    return coords

def plot2D(coords,outdir,outfile):
    
    if outdir[-1] == '/' and outfile[0] != '/' and '.' not in outfile:
        outfilepath = outdir + outfile + '.svg'
    elif outdir[-1] != '/' and outfile[0] == '/' and '.' not in outfile:
        outfilepath = outdir + outfile + '.svg'
    elif outdir[-1] != '/' and outfile[0] != '/' and '.' not in outfile:
        outfilepath = outdir + '/' + outfile + '.svg'
    elif outdir[-1] == '/' and outfile[0] == '/' and '.' not in outfile:
        outfilepath = outdir + outfile[1:] + '.svg'
    elif outdir[-1] == '/' and outfile[0] != '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + newoutfile + '.svg'
    elif outdir[-1] != '/' and outfile[0] == '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + newoutfile + '.svg'
    elif outdir[-1] != '/' and outfile[0] != '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + '/' + newoutfile + '.svg'
    elif outdir[-1] == '/' and outfile[0] == '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + newoutfile[1:] + '.svg'
    
    figure = plt.figure(1)
    for coord in coords:
        plt.plot(coord[0],coord[1],'b.')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('2D manifold')
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)
    plt.axis([1.1*mins[0], 1.1*maxs[0], 1.1*mins[1], 1.1*maxs[1]])
    plt.grid(True)
    
    figure.savefig(outfilepath)

def map2D(comp_coords,phys_coords,outdir,outfile):
    
    if outdir[-1] == '/' and outfile[0] != '/' and '.' not in outfile:
        outfilepath = outdir + outfile + '.svg'
    elif outdir[-1] != '/' and outfile[0] == '/' and '.' not in outfile:
        outfilepath = outdir + outfile + '.svg'
    elif outdir[-1] != '/' and outfile[0] != '/' and '.' not in outfile:
        outfilepath = outdir + '/' + outfile + '.svg'
    elif outdir[-1] == '/' and outfile[0] == '/' and '.' not in outfile:
        outfilepath = outdir + outfile[1:] + '.svg'
    elif outdir[-1] == '/' and outfile[0] != '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + newoutfile + '.svg'
    elif outdir[-1] != '/' and outfile[0] == '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + newoutfile + '.svg'
    elif outdir[-1] != '/' and outfile[0] != '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + '/' + newoutfile + '.svg'
    elif outdir[-1] == '/' and outfile[0] == '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + newoutfile[1:] + '.svg'
        
    figure = plt.figure(1)
    plt.title('2D map')
    sub1 = plt.subplot(1,2,1,adjustable='box',aspect=1.0)
    for coord in comp_coords:
        sub1.plot(coord[0],coord[1],'b.')
    sub1.set_xlabel(r'$\xi$')
    sub1.set_ylabel(r'$\eta$')
    sub1.set_title('Computational Domain')
    mins = comp_coords.min(axis=0)
    maxs = comp_coords.max(axis=0)
    sub1.axis([1.1*mins[0], 1.1*maxs[0], 1.1*mins[1], 1.1*maxs[1]])
    sub1.grid(True)
    sub2 = plt.subplot(1,2,2,adjustable='box',aspect=1.0)
    for coord in phys_coords:
        sub2.plot(coord[0],coord[1],'b.')
    sub2.set_xlabel(r'$x$')
    sub2.set_ylabel(r'$y$')
    sub2.set_title('Physical Domain')
    mins = phys_coords.min(axis=0)
    maxs = phys_coords.max(axis=0)
    sub2.axis([1.1*mins[0], 1.1*maxs[0], 1.1*mins[1], 1.1*maxs[1]])
    sub2.grid(True)
    
    plt.tight_layout()
    figure.savefig(outfilepath)

def plot3D(coords,outdir,outfile):
    
    if outdir[-1] == '/' and outfile[0] != '/' and '.' not in outfile:
        outfilepath = outdir + outfile + '.svg'
    elif outdir[-1] != '/' and outfile[0] == '/' and '.' not in outfile:
        outfilepath = outdir + outfile + '.svg'
    elif outdir[-1] != '/' and outfile[0] != '/' and '.' not in outfile:
        outfilepath = outdir + '/' + outfile + '.svg'
    elif outdir[-1] == '/' and outfile[0] == '/' and '.' not in outfile:
        outfilepath = outdir + outfile[1:] + '.svg'
    elif outdir[-1] == '/' and outfile[0] != '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + newoutfile + '.svg'
    elif outdir[-1] != '/' and outfile[0] == '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + newoutfile + '.svg'
    elif outdir[-1] != '/' and outfile[0] != '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + '/' + newoutfile + '.svg'
    elif outdir[-1] == '/' and outfile[0] == '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + newoutfile[1:] + '.svg'
    
    figure = plt.figure(1)
    ax = figure.add_subplot(111, projection='3d')
    ax.scatter(coords[:,0],coords[:,1],coords[:,2])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.title('3D manifold')
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)
    ax.auto_scale_xyz([1.1*mins[0], 1.1*maxs[0]], [1.1*mins[1], 1.1*maxs[1]], [1.1*mins[2], 1.1*maxs[2]])
    ax.grid(True)
    
    figure.savefig(outfilepath)

def map3D(comp_coords,phys_coords,outdir,outfile):
    
    if outdir[-1] == '/' and outfile[0] != '/' and '.' not in outfile:
        outfilepath = outdir + outfile + '.svg'
    elif outdir[-1] != '/' and outfile[0] == '/' and '.' not in outfile:
        outfilepath = outdir + outfile + '.svg'
    elif outdir[-1] != '/' and outfile[0] != '/' and '.' not in outfile:
        outfilepath = outdir + '/' + outfile + '.svg'
    elif outdir[-1] == '/' and outfile[0] == '/' and '.' not in outfile:
        outfilepath = outdir + outfile[1:] + '.svg'
    elif outdir[-1] == '/' and outfile[0] != '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + newoutfile + '.svg'
    elif outdir[-1] != '/' and outfile[0] == '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + newoutfile + '.svg'
    elif outdir[-1] != '/' and outfile[0] != '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + '/' + newoutfile + '.svg'
    elif outdir[-1] == '/' and outfile[0] == '/' and '.' in outfile:
        parts = outfile.split(".")
        newoutfile = ''.join(parts[:-1])
        outfilepath = outdir + newoutfile[1:] + '.svg'
    
    figure = plt.figure(1)
    bigAxes = plt.axes(frameon=False)     # hide frame
    plt.xticks([])                        # don't want to see any ticks on this axis
    plt.yticks([])
    cur_axes = plt.gca()
    cur_axes.axes.get_xaxis().set_visible(False)
    cur_axes.axes.get_yaxis().set_visible(False)
    plt.title('3D map')
    ax1 = figure.add_subplot(121, projection='3d',adjustable='box',aspect=1.0)
    ax1.scatter(comp_coords[:,0],comp_coords[:,1],comp_coords[:,2])
    ax1.set_xlabel(r'$\xi$')
    ax1.set_ylabel(r'$\eta$')
    ax1.set_zlabel(r'$\zeta$')
    ax1.set_title('Computational Domain')
    ax1.azim = -160
    ax1.elev = 30
    mins = comp_coords.min(axis=0)
    maxs = comp_coords.max(axis=0)
    ax1.auto_scale_xyz([1.1*mins[0], 1.1*maxs[0]], [1.1*mins[1], 1.1*maxs[1]], [1.1*mins[2], 1.1*maxs[2]])
    ax1.grid(True)
    
    ax2 = figure.add_subplot(122, projection='3d',adjustable='box',aspect=1.0)
    ax2.scatter(phys_coords[:,0],phys_coords[:,1],phys_coords[:,2])
    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$y$')
    ax2.set_zlabel(r'$z$')
    ax2.set_title('Physical Domain')
    ax2.azim = -160
    ax2.elev = 30
    mins = phys_coords.min(axis=0)
    maxs = phys_coords.max(axis=0)
    ax2.auto_scale_xyz([1.1*mins[0], 1.1*maxs[0]], [1.1*mins[1], 1.1*maxs[1]], [1.1*mins[2], 1.1*maxs[2]])
    ax2.grid(True)
    
    plt.tight_layout()
    figure.savefig(outfilepath)

def main(argv):
    
    # Read the command line, throw error if not option is provided
    try:
        opts, args = getopt.getopt(argv,'hD:x:y:z:t:o:d:c:w:',['help','Help','dim','dimension','Dimension','Dim','Nx','Ny','Nz','test','Test','out','output','outfile','dir','directory','outdir','outdirectory','cppfile','cppdir'])
    except getopt.GetoptError:
        print 'generate_manifolds.py -D <dimension> -x <number of points in x-direction> -y <number of points in y-direction> -z <number of points in z-direction> -t <test case>'
        print '                      -o <output file> -d <output directory> -c <c++ output file> -w <c++ work directory>'
        sys.exit(2)
    
    # Parse the options and create corresponding variables
    for opt, arg in opts:
        if opt in ('-h', '--help','--Help'):
            print ' '  
            print ' '  
            print '********************************************************************************* '  
            print ' '  
            print '                                     EUSMAT'  
            print '                         European School of Materials'  
            print ' '  
            print '                                     DocMASE'  
            print '                  DOCTORATE IN MATERIALS SCIENCE AND ENGINEERING'  
            print ' '  
            print ' '  
            print '       MECHANICS OF EXTREME THIN COMPOSITE LAYERS FOR AEROSPACE APPLICATIONS'  
            print ' '  
            print ' '  
            print '                         Differential geometry tools'  
            print ' '
            print '                          Generation of test cases'
            print ' '
            print '                                      by'  
            print ' '  
            print '                             Luca Di Stasio, 2016'  
            print ' '  
            print ' '  
            print '********************************************************************************* '  
            print ''
            print 'Program syntax:'
            print 'generate_manifolds.py -D <dimension> -x <number of points in x-direction> -y <number of points in y-direction> -z <number of points in z-direction> -t <test case>'
            print '                      -o <output file> -d <output directory> -c <c++ output file> -w <c++ work directory>'
            print ''
            print 'Mandatory arguments:'
            print '-x <number of points in x-direction> -y <number of points in y-direction> -t <test case>'
            print ''
            print 'Optional arguments:'
            print '-D <dimension> -z <number of points in z-direction> -o <output file (without extension)> -d <output directory> -c <c++ output file> -w <c++ work directory>'
            print ''
            print 'Default values:'
            print '-D <dimension>                       ===> 2'
            print '-o <output file (without extension)> ===> example'
            print '-d <output directory>                ===> /home/ubuntu/workspace/Differential_geometry/WD/Input_files'
            print '-c <c++ output file>                 ===> geometry_benchmark'
            print '-w <c++ work directory>              ===> Output_files'
            print ''
            print 'Available test cases:'
            print ''
            print '-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
            print '| Case                                           Command-line option                                                                                                                                            |'
            print '-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
            print ''
            print '2D planar, equal squares                         c1d2l<side length>cx<x coordinate of center>cy<y coordinate of center>'
            print 'Case 1 (c1). Two-dimensional (d2) test case. Two equal squares of length l and centered in (cx,cy) mapped unto each other.'
            print ''
            print '2D planar, equal rectangles                      c2d2lx<x side length>ly<y side length>cx<x coordinate of center>cy<y coordinate of center>'
            print 'Case 2 (c2). Two-dimensional (d2) test case. Two equal rectangles with sides of length lx and ly (respectively, the side parallel to the x-axis and the one parallel to the y-axis) and centered in (cx,cy) mapped unto each other.'
            print ''
            print '2D planar, stretched square                      c3d2l<side length>cx<x coordinate of center>cy<y coordinate of center>s<stretch factor>'
            print 'Case 3 (c3). Two-dimensional (d2) test case. A square of length l and centered in (cx,cy) mapped unto its image isotropically stretched by a factor of s.'
            print ''
            print '2D planar, stretched rectangle                   c4d2lx<x side length>ly<y side length>cx<x coordinate of center>cy<y coordinate of center>s<stretch factor>'
            print 'Case 4 (c4). Two-dimensional (d2) test case. A rectangle with sides of length lx and ly (respectively, the side parallel to the x-axis and the one parallel to the y-axis) and centered in (cx,cy)'
            print 'mapped unto its image isotropically stretched by a factor of s.'
            print ''
            print '2D planar, rounded square                        c5d2l<side length>cx<x coordinate of center>cy<y coordinate of center>r<rounding radius>crx<x coordinate of rounding center>cry<y coordinate of rounding center>'
            print 'Case 5 (c5). Two-dimensional (d2) test case. A square of length l and centered in (cx,cy) mapped unto its image with sides replaced by circular arcs of center (crx, cry) with respect to the square center'
            print '(thus absolute coordinates equal to cx+crx and cy+cry) and radius r.'
            print ''
            print '2D planar, rounded rectangle                     c6d2lx<x side length>ly<y side length>cx<x coordinate of center>cy<y coordinate of center>rx<rounding radius for x side>ry<rounding radius for y side>'
            print '                                                 crxx<x coordinate of rounding center for x side>crxy<y coordinate of rounding center for x side>'
            print '                                                 cryx<x coordinate of rounding center for y side>cryy<y coordinate of rounding center for y side>'
            print 'Case 6 (c6). Two-dimensional (d2) test case.  A rectangle with sides of length lx and ly (respectively, the side parallel to the x-axis and the one parallel to the y-axis) and centered in (cx,cy)'
            print 'mapped unto its image with sides replaced by circular arcs of center (crxx, crxy) and (cryx, cryy) (respectively for the side parallel to the x-axis and the one parallel to the y-axis)'
            print 'with respect to the rectangle center (thus absolute coordinates equal to -cx+crxx,cy+crxy- and -cx+cryx,cy+cryy-) and radiuses rx and ry.'
            print ''
            print '2D planar, circle                                c7d2l<side length>cx<x coordinate of center>cy<y coordinate of center>r<radius>'
            print 'Case 7 (c7). Two-dimensional (d2) test case. A square of length l and centered in (cx,cy) mapped unto a circle of radius r with the same center.'
            print ''
            print '2D planar, ellipse                               c8d2l<side length>cx<x coordinate of center>cy<y coordinate of center>a<major semi-axis>b<minor semi-axis>'
            print 'Case 8 (c8). Two-dimensional (d2) test case. A square of length l and centered in (cx,cy) mapped unto an ellipse with major and minor semi-axes respectively equal to a and b with the same center.'
            print ''
            print '3D, equal cubes                                  c9d3l<side length>cx<x coordinate of center>cy<y coordinate of center>cz<z coordinate of center>'
            print 'Case 9 (c9). Three-dimensional (d3) test case. Two equal cubes of length l and centered in (cx,cy,cz) mapped unto each other.'
            print ''
            print '3D, equal parallelepipeds                        c10d3lx<x side length>ly<y side length>lz<z side length>cx<x coordinate of center>cy<y coordinate of center>cz<z coordinate of center>'
            print 'Case 10 (c10). Three-dimensional (d3) test case. Two equal parallelepipeds with sides of length lx, ly and lz (respectively, the side parallel to the x-axis, y-axis and z-axis) and centered in (cx,cy,cz) mapped unto each other.'
            print ''
            print '3D, stretched cube                               c11d3l<side length>cx<x coordinate of center>cy<y coordinate of center>cz<z coordinate of center>s<stretch factor>'
            print 'Case 11 (c11). Three-dimensional (d3) test case. A cube of length l and centered in (cx,cy,cz) mapped unto its image isotropically stretched by a factor of s.'
            print ''
            print '3D, stretched parallelepiped                     c12d3lx<x side length>ly<y side length>lz<z side length>cx<x coordinate of center>cy<y coordinate of center>cz<z coordinate of center>s<stretch factor>'
            print 'Case 12 (c12). Three-dimensional (d3) test case. A parallelepiped of length lx, ly and lz (respectively, the side parallel to the x-axis, y-axis and z-axis) and centered in (cx,cy,cz) mapped unto its image isotropically stretched by a factor of s.'
            print ''
            print '3D, pipe with rounded square section             c13d3lx<x side length>ly<y and z side length>cx<x coordinate of center>cy<y coordinate of center>cz<z coordinate of center>r<rounding radius>crx<x coordinate of rounding center>cry<y coordinate of rounding center>'
            print 'Case 13 (c13). Three-dimensional (d3) test case. A parallelepiped of length lx, ly and again ly (respectively, the side parallel to the x-axis, y-axis and z-axis) and centered in (cx,cy,cz) mapped unto a pipe with rounded square section'
            print 'with sides replaced by circular arcs of center (cry, crz) with respect to the square center (thus absolute coordinates equal to cy+cry and cz+crz) and radius r. The pipe axis is in the x-direction.'
            print ''
            print '3D, cylinder                                     c14d3lx<x side length>ly<y and z side length>cx<x coordinate of center>cy<y coordinate of center>cz<z coordinate of center>r<cylinder radius>'
            print 'Case 14 (c14). Three-dimensional (d3) test case. A parallelepiped of length lx, ly and again ly (respectively, the side parallel to the x-axis, y-axis and z-axis) and centered in (cx,cy,cz) mapped unto a cylinder of radius r'
            print 'with the same center. The cylinder axis is in the x-direction. '
            print ''
            print '3D, sphere                                       c15d3l<side length>cx<x coordinate of center>cy<y coordinate of center>cz<z coordinate of center>r<radius>'
            print 'Case 15 (c15). Three-dimensional (d3) test case. A cube of length l and centered in (cx,cy,cz) mapped unto unto a sphere of radius r with the same center.'
            print ''
            print '3D, ellipsoid                                    c16d3l<side length>cx<x coordinate of center>cy<y coordinate of center>cz<z coordinate of center>a<first principal semi-axis>b<first principal semi-axis>c<first principal semi-axis>'
            print 'Case 16 (c16). Three-dimensional (d3) test case. A cube of length l and centered in (cx,cy,cz) mapped unto an ellipsoid with principal semi-axes respectively equal to a, b and c with the same center.'
            print ''
            sys.exit()
        elif opt in ("-D", "--dim","--dimension","--Dimension","--Dim"):
            dim = int(arg)
        elif opt in ("-x", "--Nx"):
            Nx = int(arg)
        elif opt in ("-y", "--Ny"):
            Ny = int(arg)
        elif opt in ("-z", "--Nz"):
            Nz = int(arg)
        elif opt in ("-t", '--test','--Test'):
            test = arg
        elif opt in ("-o", '--output','--outfile'):
            outfilename = arg
        elif opt in ("-d", '--dir','--directory','--outdir','--outdirectory'):
            outdir = arg
        elif opt in ("-c", '--cppfile'):
            cppfile = arg
        elif opt in ("-w", '--cppdir'):
            cppdir = arg
    
    # Check the existence of variables: if a required variable is missing, an error is thrown and program is terminated; if an optional variable is missing, it is set to the default value
    if 'dim' not in locals():
        dim = 2
    if 'Nx' not in locals():
        print 'Error: option x not provided.'
        sys.exit(2)
    if 'Ny' not in locals():
        print 'Error: option y not provided.'
        sys.exit(2)
    if 'Nz' not in locals() and dim == 3:
        print 'Error: option z not provided.'
        sys.exit(2)
    elif 'Nz' not in locals() and dim == 2:
        Nz = 0
    if 'test' not in locals():
        print 'Error: option test not provided.'
        sys.exit(2)
    if 'outfilename' not in locals():
        outfilename = 'example'
    if 'outdir' not in locals():
        outdir = '/home/ubuntu/workspace/Differential_geometry/WD/Input_files'
    if 'cppfile' not in locals():
        cppfile = 'geometry_benchmark'
    if 'cppdir' not in locals():
        cppdir = 'Output_files'
        
    # Parse the test option and assign values to corresponding parameters
    descriptors = re.split('(\d+)', test)
    for index, descriptor in enumerate(descriptors):
        if descriptor == 'c':
            case = int(descriptors[index+1])
        elif descriptor == 'l':
            l = float(descriptors[index+1])
        elif descriptor == 'lx':
            lx = float(descriptors[index+1])
        elif descriptor == 'ly':
            ly = float(descriptors[index+1])
        elif descriptor == 'lz':
            lz = float(descriptors[index+1])
        elif descriptor == 'cx':
            cx = float(descriptors[index+1])
        elif descriptor == 'cy':
            cy = float(descriptors[index+1])
        elif descriptor == 'cz':
            cz = float(descriptors[index+1])
        elif descriptor == 'crx':
            crx = float(descriptors[index+1])
        elif descriptor == 'cry':
            cry = float(descriptors[index+1])
        elif descriptor == 'crz':
            crz = float(descriptors[index+1])
        elif descriptor == 'crxx':
            crxx = float(descriptors[index+1])
        elif descriptor == 'crxy':
            crxy = float(descriptors[index+1])
        elif descriptor == 'cryx':
            cryx = float(descriptors[index+1])
        elif descriptor == 'cryy':
            cryy = float(descriptors[index+1])
        elif descriptor == 'r':
            r = float(descriptors[index+1])
        elif descriptor == 'rx':
            rx = float(descriptors[index+1])
        elif descriptor == 'ry':
            ry = float(descriptors[index+1])
        elif descriptor == 's':
            s = float(descriptors[index+1])
        elif descriptor == 'a':
            a = float(descriptors[index+1])
        elif descriptor == 'b':
            b = float(descriptors[index+1])
        elif descriptor == 'c':
            c = float(descriptors[index+1])
    
    if 'case' not in locals():
        print 'Error: option provided to -t <test> was wrongly formatted. Case not specified.'
        sys.exit(2)

    if case == 1:
        if 'l' not in locals() or 'cx' not in locals() or 'cy' not in locals() or 'Nx' not in locals() or 'Ny' not in locals():
            print 'Error: a geometric parameter was not provided. Check the command line option you submitted (-t option) for errors.'
            sys.exit(2)
        else:
            comp_coords = square(l,cx,cy,Nx,Ny)
            phys_coords = square(l,cx,cy,Nx,Ny)
            map2D(comp_coords,phys_coords,outdir,outfilename)
            write_csv(dim,comp_coords,outdir,outfilename+'comp')
            write_csv(dim,phys_coords,outdir,outfilename+'phys')
            write_inp([dim,Nx,Ny,Nz,outfilename+'comp',outfilename+'phys',cppdir,cppfile],outdir,outfilename)
    elif case == 2:
        if 'lx' not in locals() or 'ly' not in locals() or 'cx' not in locals() or 'cy' not in locals() or 'Nx' not in locals() or 'Ny' not in locals():
            print 'Error: a geometric parameter was not provided. Check the command line option you submitted (-t option) for errors.'
            sys.exit(2)
        else:
            comp_coords = rectangle(lx,ly,cx,cy,Nx,Ny)
            phys_coords = rectangle(lx,ly,cx,cy,Nx,Ny)
            map2D(comp_coords,phys_coords,outdir,outfilename)
            write_csv(dim,comp_coords,outdir,outfilename+'comp')
            write_csv(dim,phys_coords,outdir,outfilename+'phys')
            write_inp([dim,Nx,Ny,Nz,outfilename+'comp',outfilename+'phys',cppdir,cppfile],outdir,outfilename)
    elif case == 3:
        if 'l' not in locals() or 'cx' not in locals() or 'cy' not in locals() or 'Nx' not in locals() or 'Ny' not in locals() or 's' not in locals():
            print 'Error: a geometric parameter was not provided. Check the command line option you submitted (-t option) for errors.'
            sys.exit(2)
        else:
            comp_coords = square(l,cx,cy,Nx,Ny)
            phys_coords = square(s*l,cx,cy,Nx,Ny)
            map2D(comp_coords,phys_coords,outdir,outfilename)
            write_csv(dim,comp_coords,outdir,outfilename+'comp')
            write_csv(dim,phys_coords,outdir,outfilename+'phys')
            write_inp([dim,Nx,Ny,Nz,outfilename+'comp',outfilename+'phys',cppdir,cppfile],outdir,outfilename)
    elif case == 4:
        if 'lx' not in locals() or 'ly' not in locals() or 'cx' not in locals() or 'cy' not in locals() or 'Nx' not in locals() or 'Ny' not in locals() or 's' not in locals():
            print 'Error: a geometric parameter was not provided. Check the command line option you submitted (-t option) for errors.'
            sys.exit(2)
        else:
            comp_coords = rectangle(lx,ly,cx,cy,Nx,Ny)
            phys_coords = rectangle(s*lx,s*ly,cx,cy,Nx,Ny)
            map2D(comp_coords,phys_coords,outdir,outfilename)
            write_csv(dim,comp_coords,outdir,outfilename+'comp')
            write_csv(dim,phys_coords,outdir,outfilename+'phys')
            write_inp([dim,Nx,Ny,Nz,outfilename+'comp',outfilename+'phys',cppdir,cppfile],outdir,outfilename)
    elif case == 5:
        print 'Feature still not available, sorry!'
    elif case == 6:
        print 'Feature still not available, sorry!'
    elif case == 7:
        print 'Feature still not available, sorry!'
    elif case == 8:
        print 'Feature still not available, sorry!'
    elif case == 9:
        if 'l' not in locals() or 'cx' not in locals() or 'cy' not in locals() or 'cz' not in locals() or 'Nx' not in locals() or 'Ny' not in locals() or 'Nz' not in locals():
            print 'Error: a geometric parameter was not provided. Check the command line option you submitted (-t option) for errors.'
            sys.exit(2)
        else:
            comp_coords = cube(l,cx,cy,cz,Nx,Ny,Nz)
            phys_coords = cube(l,cx,cy,cz,Nx,Ny,Nz)
            map3D(comp_coords,phys_coords,outdir,outfilename)
            write_csv(dim,comp_coords,outdir,outfilename+'-comp')
            write_csv(dim,phys_coords,outdir,outfilename+'-phys')
            write_inp([dim,Nx,Ny,Nz,outfilename+'-comp',outfilename+'-phys',cppdir,cppfile],outdir,outfilename)
    elif case == 10:
        if 'lx' not in locals() or 'ly' not in locals() or 'lz' not in locals() or 'cx' not in locals() or 'cy' not in locals() or 'cz' not in locals() or 'Nx' not in locals() or 'Ny' not in locals() or 'Nz' not in locals():
            print 'Error: a geometric parameter was not provided. Check the command line option you submitted (-t option) for errors.'
            sys.exit(2)
        else:
            comp_coords = parallelepiped(lx,ly,lz,cx,cy,cz,Nx,Ny,Nz)
            phys_coords = parallelepiped(lx,ly,lz,cx,cy,cz,Nx,Ny,Nz)
            map3D(comp_coords,phys_coords,outdir,outfilename)
            write_csv(dim,comp_coords,outdir,outfilename+'-comp')
            write_csv(dim,phys_coords,outdir,outfilename+'-phys')
            write_inp([dim,Nx,Ny,Nz,outfilename+'-comp',outfilename+'-phys',cppdir,cppfile],outdir,outfilename)
    elif case == 11:
        if 'l' not in locals() or 'cx' not in locals() or 'cy' not in locals() or 'cz' not in locals() or 'Nx' not in locals() or 'Ny' not in locals() or 'Nz' not in locals() or 's' not in locals():
            print 'Error: a geometric parameter was not provided. Check the command line option you submitted (-t option) for errors.'
            sys.exit(2)
        else:
            comp_coords = cube(l,cx,cy,cz,Nx,Ny,Nz)
            phys_coords = cube(s*l,cx,cy,cz,Nx,Ny,Nz)
            map3D(comp_coords,phys_coords,outdir,outfilename)
            write_csv(dim,comp_coords,outdir,outfilename+'-comp')
            write_csv(dim,phys_coords,outdir,outfilename+'-phys')
            write_inp([dim,Nx,Ny,Nz,outfilename+'-comp',outfilename+'-phys',cppdir,cppfile],outdir,outfilename)
    elif case == 12:
        if 'lx' not in locals() or 'ly' not in locals() or 'lz' not in locals() or 'cx' not in locals() or 'cy' not in locals() or 'cz' not in locals() or 'Nx' not in locals() or 'Ny' not in locals() or 'Nz' not in locals() or 's' not in locals():
            print 'Error: a geometric parameter was not provided. Check the command line option you submitted (-t option) for errors.'
            sys.exit(2)
        else:
            comp_coords = parallelepiped(lx,ly,lz,cx,cy,cz,Nx,Ny,Nz)
            phys_coords = parallelepiped(s*lx,s*ly,s*lz,cx,cy,cz,Nx,Ny,Nz)
            map3D(comp_coords,phys_coords,outdir,outfilename)
            write_csv(dim,comp_coords,outdir,outfilename+'-comp')
            write_csv(dim,phys_coords,outdir,outfilename+'-phys')
            write_inp([dim,Nx,Ny,Nz,outfilename+'-comp',outfilename+'-phys',cppdir,cppfile],outdir,outfilename)
    elif case == 13:
        print 'Feature still not available, sorry!'
    elif case == 14:
        print 'Feature still not available, sorry!'
    elif case == 15:
        print 'Feature still not available, sorry!'
    elif case == 16:
        print 'Feature still not available, sorry!'    
    


if __name__ == "__main__":
    main(sys.argv[1:])