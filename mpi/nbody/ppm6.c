/*******************************************************************************
 *
 * ppm6.c
 * by Ben Witzen <b.a.witzen@uva.nl>
 * November 12th, 2015
 *
 * The image I/O functions as provided by the staff.
 *
 ******************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include "ppm6.h"

unsigned char *map_P6(char *filename, int *xdim, int *ydim)
{
  int	fsize;
  unsigned char	*map;

	/* The following is a fast and sloppy way to
	   map a color raw PPM (P6) image file
	*/
	int fd;
	unsigned char *p;
	int maxval;

	/* First, open the file... */
	if ((fd = open(filename, O_RDWR)) < 0) {
		return((unsigned char *) 0);
	}

	/* Read size and map the whole file... */
	fsize = lseek(fd, ((off_t) 0), SEEK_END);
	map = ((unsigned char *)
	       mmap(0,		/* Put it anywhere */
		    fsize,	/* Map the whole file */
		    (PROT_READ|PROT_WRITE),	/* Read/write */
		    MAP_SHARED,	/* Not just for me */
		    fd,		/* The file */
		    0));	/* Right from the start */
	if (map == ((unsigned char *) -1)) {
		close(fd);
		return((unsigned char *) 0);
	}

	/* File should now be mapped; read magic value */
	p = map;
	if (*(p++) != 'P') goto ppm_exit;
	switch (*(p++)) {
	case '6':
		break;
	default:
		goto ppm_exit;
	}

#define	Eat_Space \
	while ((*p == ' ') || \
	       (*p == '\t') || \
	       (*p == '\n') || \
	       (*p == '\r') || \
	       (*p == '#')) { \
		if (*p == '#') while (*(++p) != '\n') ; \
		++p; \
	}

	Eat_Space;		/* Eat white space and comments */

#define	Get_Number(n) \
	{ \
		int charval = *p; \
 \
		if ((charval < '0') || (charval > '9')) goto ppm_exit; \
 \
		n = (charval - '0'); \
		charval = *(++p); \
		while ((charval >= '0') && (charval <= '9')) { \
			n *= 10; \
			n += (charval - '0'); \
			charval = *(++p); \
		} \
	}

	Get_Number(*xdim);	/* Get image width */

	Eat_Space;		/* Eat white space and comments */
	Get_Number(*ydim);	/* Get image width */

	Eat_Space;		/* Eat white space and comments */
	Get_Number(maxval);	/* Get image max value */

	/* Should be 8-bit binary after one whitespace char... */
	if (maxval > 255) {
ppm_exit:
		close(fd);
		munmap(map, fsize);
		return((unsigned char *) 0);
	}
	if ((*p != ' ') &&
	    (*p != '\t') &&
	    (*p != '\n') &&
	    (*p != '\r')) goto ppm_exit;

	/* Here we are... next byte begins the 24-bit data */
	return(p + 1);

	/* Notice that we never clean-up after this:

	   close(fd);
	   munmap(map, fsize);

	   However, this is relatively harmless;
	   they will go away when this process dies.
	*/
}

#undef	Eat_Space
#undef	Get_Number
