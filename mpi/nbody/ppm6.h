/*******************************************************************************
 *
 * ppm6.h
 * by Ben Witzen <b.a.witzen@uva.nl>
 * November 12th, 2015
 *
 * The image I/O functions as provided by the staff.
 *
 ******************************************************************************/

#ifndef PPM6_H
#define PPM6_H

#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>

unsigned char *map_P6(char *filename, int *xdim, int *ydim);

#endif
