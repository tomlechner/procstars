// Star Catalog Image Generator
// By Tom Lechner, tomlechner.com
// 2013
//
// This program creates star globes from the wonderful Tycho2 star catalog,
// and also the Principal Galaxy Catalog.
//
// Usage:
//  ./procstars  -c ../Tycho.dat  -w32768  -H 7 -m12 -S 10 -a 40 -g
//  ./procstars  -c ../Tycho.cat  -S 6 -m 12 -s 12 -w 32768
//  ./procstars  -c ../Tycho.cat  -S 6 -m 12 -s 6 -w 16384
//  ./procstars  -c ../Tycho.cat  -o stars.tif
//
//
//
// This is released with the MIT license.
// This was first inspired by a generator written by Nathan Bergey from 01-22-2010,
// which he released CC Attribution 3.0 as Processing code, which processed the 
// Tycho 2 catalog.
//
// Fyi, the ISS takes up about 30 arc seconds.
// The moon takes up about 30 arc minutes.
// at 8192 px wide, 23 pixels make 1 degree. 1/3 of a pixel is one arc minute, 1/200 of a pixel is 1 arc second.
// at 32768,        91 pixels make 1 degree, 1.5 pixels make one arc minute. 
//
//
// The MIT License (MIT)
// Copyright (c) 2013 Tom Lechner
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
//
//


//Todo:
//  implement galaxy approximations with principal galaxy catalog's morphology field
//  project billboard stuff like nebulae, LMC and SMC, M31, space stations, etc
//  make halo scale with magnitude
//  find nebula catalog?
//  what are those huge red smears in the sky??
//  fake the big galaxies LMC+SMC, m31?



//Tycho information:
//
//Dimmest magnitude:   0
//Brightest magnitude: 15
//Magnitudes:
//  0: 48
//  1: 52
//  2: 252
//  3: 677
//  4: 1850
//  5: 4774
//  6: 12665
//  7: 30248
//  8: 79522
//  9: 207629
//  10: 533939
//  11: 1121779
//  12: 533238
//  13: 13118
//  14: 121
//  15: 1
//
//For the tycho catalog, in principal you can see all the stars in the catalog,
//but after magnitude of about 12.1, you begin to see what appear to be striations 
//from where the satellite did not scan.
//



#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <stdint.h>
#include <getopt.h>

#include </usr/include/GraphicsMagick/Magick++.h>

//#include <lax/laxoptions.h>



#define PI (M_PI)
#define TWO_PI (2*M_PI)

//map: scale from [omin,omax] to [nmin,nmax]
#define map(value, omin,omax, nmin,nmax) (nmin+(value-omin)/(omax-omin)*(nmax-nmin))
#define constrain(v, min,max) (v<min ? min : (v>max ? max : v))
#define radians(deg) (deg/180.*M_PI)
#define degrees(rad) (rad*180./M_PI)



using namespace std;
using namespace Magick;


//options:
const char *filename="stars.tif";
long width=8192;
long height=4096;
double magnitude=15.1;
double maxmagnitude=-5;
int galactic=0;
double maxstarsize=3; //pixels of biggest star
double bigthreshhold=4; //magnitude < this get big treatment
int transparent=0;
double alphaamp=0;
double usehalo=0;


//stuff to make computations easier:
double pixelwidth=1./width;
unsigned char *data=NULL;
unsigned char *halo=NULL;
int halowidth=0;


//----------------------- Encapsulate a catalog (work in progress) -------------------------
class RenderContext
{
  public:
	double min_asc;
	double max_asc;
	double min_dec;
	double max_dec;

	long width;
	long height;
	unsigned char *data;
};

enum CatalogTypes
{
	Tycho2,
	PGC
};

class Catalog
{
  public:
	char *name;
	char *filename;
	CatalogTypes type; //Tycho2, PGC, point cloud (generated far off galaxies), sphere cloud (surrounding earth)
	
	int *magnitude_distribution;
	int nummags;
	double minimum_magnitude;
	double maximum_magnitude;

	int min_mag_cutoff;
	int max_mag_cutoff;

	Catalog(const char *nname, const char *nfile, CatalogTypes ntype);
	virtual ~Catalog();
	virtual int Render(RenderContext *context);
	virtual int GetStats();
};





//------------------------- forward decs -----------------------------
class flatvector
{
  public:
	double x,y;
	flatvector(double xx,double yy) { x=xx; y=yy; }
};
flatvector Eq2Gal(float ra, float dec);
void indexToRgb(float index, float vmag, ColorRGB &color);
double indexToRed(float index);
double indexToGreen(float index);
double indexToBlue(float index);
void drawStar(float ra, float dec, float vmag, float bmag);
void CreateStockHalo(int w,double halosize, unsigned char *halo);
double dms(const char *pos);
double hms(const char *pos);





void help()
{
	cout << endl;
	cout << "Star Catalog Image Generator"<<endl;
	cout << "By Tom Lechner"<<endl;
	cout << "Version 0.000000001"<<endl;
	cout << "\nprocstars [options]\n\n";

	cout << "--width, -w 8192            Width of equirectangular texture image"<<endl;
	cout << "--output-file, -o file      Filename for generated texture image (will save as tif)"<<endl;
	cout << "--catalog-file, -c catalog  File containing the Tycho 2 Star Catalog"<<endl;
	cout << "--pgc-file, -p catalog      File containing the Principal Galaxy Catalog"<<endl;
	cout << "--magnitude, -m 15.1        Minimum brightness. stars in catalog are 0 (bright)\n"
			"                              to 15 (very dim). 6 is dimmest to human naked eye."<<endl;
	cout << "--transparent, -T           Put stars on transparency rather than a black background"<<endl;
	cout << "--maxstarsize, -s 3         Pixels wide of brightest star"<<endl;
	cout << "--maxstarangle, -a 3        Arc minutes wide of brightest star"<<endl;
	cout << "--bigthreshhold, -S 4       Magnitude brighter than this gets special size treatment"<<endl;
	cout << "--alpha-amp, -A 0           Adjustment to artifically pump up dim stars. Between 0 and 1"<<endl;
	cout << "--galactic, -g              Make it so the Milky Way is horizontal"<<endl;
	cout << "--halo, -H 1.5              Render halos on bright stars, number is ratio of \n"
		    "                              halo diameter to star diameter"<<endl;
	exit(0);
}

int main(int argc, char **argv)
{
	InitializeMagick(*argv);

	const char *tycho_file=NULL;
	//const char *tycho_file="Tycho.dat";
	const char *pgc_file=NULL;
	//const char *pgc_file="Tycho.dat";

     // parse options
    static struct option long_options[] = {
            { "width",          1, 0, 'w' },
            { "output-file",    1, 0, 'o' },
            { "catalog-file",   1, 0, 'c' },
            { "catalog-file",   1, 0, 'p' },
            { "magnitude",      1, 0, 'm' },
            { "maxstarsize",    1, 0, 's' },
            { "maxstarangle",   1, 0, 'a' },
            { "bigthreshhold",  1, 0, 'S' },
            { "galactic",       0, 0, 'g' },
            { "halo",           1, 0, 'H' },
            { "transparent",    0, 0, 'T' },
            { "alpha-amp",      1, 0, 'A' },
			{ "version",        0, 0, 'v' },
            { "help",           0, 0, 'h' },
            { 0,0,0,0 }
        };
    int cc,index;

	double ang=-1;
    while (1) {
        cc=getopt_long(argc,argv,":w:o:c:p:m:s:S:a:A:H:gTvh",long_options,&index);
        if (cc==-1) break;
        switch(cc) {
            case ':': cerr <<"Missing parameter..."<<endl; exit(1); // missing parameter
            case '?': cerr <<"Unknown option"<<endl; exit(1);  // unknown option
            case 'h': help();    // Show usage summary, then exit
            case 'v': help(); // Show version info, then exit

            case 'w': {
				width=strtol(optarg,NULL,10);
				if (width<2) { cerr <<"Badth width value!"<<endl; exit(1); }
				height=width/2;
              } break;

            case 'o': {
				filename=optarg;
              } break;

            case 'p': {
				pgc_file=optarg;
              } break;

            case 'c': {
				tycho_file=optarg;
              } break;

            case 'm': {
				magnitude=strtod(optarg,NULL);
              } break;

            case 'g': {
				galactic=1;
              } break;

            case 'T': {
				transparent=1;
              } break;

            case 'H': {
				usehalo=strtod(optarg,NULL);
              } break;

            case 's': {
				maxstarsize=strtod(optarg,NULL);
              } break;

            case 'a': {
				ang=strtod(optarg,NULL);
              } break;

            case 'S': {
				bigthreshhold=strtod(optarg,NULL);
              } break;

            case 'A': {
				alphaamp=strtod(optarg,NULL);
              } break;

        }
    }
	if (ang>0) {
		maxstarsize=width*ang/60./360.;
	}
	if (usehalo>1) {
		halowidth=maxstarsize*usehalo;
		halo=new unsigned char[halowidth*halowidth];
		CreateStockHalo(halowidth, usehalo, halo);
	}




	pixelwidth=1./width;


	cout <<endl;
	cout << "      Tycho catalog: "<<(tycho_file?tycho_file:"none")<<endl;
	cout << "     Galaxy catalog: "<<(pgc_file?pgc_file:"none")<<endl;
	cout << "      Outputting to: "<<filename<<endl;
	cout << "         Dimensions: "<<width<<" x "<<height<<endl;
	cout << "Brightness at least: "<<magnitude<<endl;
	cout << "        Halos after: "<<bigthreshhold<<endl;





	 //allocate star image data
	unsigned char *ddata=new unsigned char[(long)width*height*4];

	int stride=width*4;
	data=ddata;
	memset(data,0,width*height*4);
	if (!transparent) {
		for (int y=0; y<height; y++) {
		  for (int x=0; x<width; x++) {
			//data[y*stride+x*4+0]=255;
			//data[y*stride+x*4+1]=255;
			//data[y*stride+x*4+2]=255;
			data[y*stride+x*4+3]=255;//opaque
		  }
		}
	}





    double Ra;   //right ascension (longitude)
    double Dec;  //declination (latitude)
    double vmag; //visual mag
    double bmag; //blue mag

	//double colorindex;
	//double colorindex_min= 10000;
	//double colorindex_max=-10000;

	int mags[50];
	int magmin=1000, magmax=-1000;
	memset(mags,0,50*sizeof(int));

	int numstars=0;




	 //--------------- Principal Galaxy Catalog processing -------------------
	int numgalaxies=0;
	if (pgc_file) {
		char *line=NULL;
		size_t n=0;
		ssize_t c;
		FILE *f=fopen(pgc_file,"r");
		if (!f) {
			cerr <<" --Fail!-- Could not open Principal Galaxy Catalog: "<<pgc_file<<endl;
			exit(1);
		}
		do {

			c=getline(&line,&n,f);
			if (c<0) break;
			if (c==0) continue;
			if (line[6]== ' ') continue; //no coordinate data on this line
			if (line[59]==' ') continue; //no magnitude data on this line

			 //the catalog lines are delimited by '|' characters, but also arranged on strict byte widths
			Ra   = hms(line+6);          // RA hours-min-sec,        bytes 7-14
			Dec  = dms(line+14);         //DEC degrees-min-sec,      bytes 15-21
			vmag = strtod(line+59,NULL); //overall visual magnitude, bytes 60-63
			bmag = vmag+.2; //makes index refer to a slightly blue object

			// *** type of galaxy, bytes 40-43
			// *** major axis, arcmin, bytes 44-49
			// *** minor axis, arcmin, bytes 52-56
			// *** Position Angle from North eastward, bytes 74-76
		 
			if (vmag>magnitude || vmag<maxmagnitude) continue;
			if (Ra==0 && Dec==0) continue;

			if (vmag<magmin) magmin=vmag;
			if (vmag>magmax) magmax=vmag;
			mags[int(vmag-magmin)]++;


			drawStar(Ra, Dec, vmag, bmag);

			numgalaxies++;
			if (numgalaxies%10000==0) cout <<"+\n";

		} while (!feof(f));
		if (line) free(line);

		fclose(f);
	}



	 //--------------- Tycho star catalog processing ---------------
	if (tycho_file) {
		char *line=NULL;
		size_t n=0;
		ssize_t c;
		FILE *f=fopen(tycho_file,"r");
		if (!f) {
			cerr <<" --Fail!-- Could not open Tycho catalog: "<<tycho_file<<endl;
			exit(1);
		}
		do {

			c=getline(&line,&n,f);
			if (c<0) break;
			if (c==0) continue;


			 //the catalog lines are delimited by '|' characters, but also arranged on strict byte widths
			Ra   = strtod(line+15,NULL);  //pieces[2],  mRA degrees,   bytes 16-27
			Dec  = strtod(line+28,NULL);  //pieces[3],  mDE degrees,   bytes 29-40
			bmag = strtod(line+110,NULL); //pieces[17], BT, tycho 2 blue magnitude,   bytes 111-116
			vmag = strtod(line+123,NULL); //pieces[19], VT, tycho 2 visual magnitude, bytes 124-129
		 
			if (vmag>magnitude || vmag<maxmagnitude) continue;
			if (Ra==0 && Dec==0) continue;

			if (vmag<magmin) magmin=vmag;
			if (vmag>magmax) magmax=vmag;
			mags[int(vmag-magmin)]++;

			 //find stats about colors:
			//colorindex=bmag-vmag;
			//if (colorindex<colorindex_min) colorindex_min=colorindex;
			//if (colorindex>colorindex_max) colorindex_max=colorindex;


			drawStar(Ra, Dec, vmag, bmag);

			numstars++;
			if (numstars%100000==0) cout <<".\n";

		} while (!feof(f));
		if (line) free(line);

		fclose(f);
	}






	//-------------- All done! print summary------------------

	 //statistics
	cout <<endl;
	cout <<"Magnitudes:"<<endl;
	for (int c=magmin; c<=magmax; c++) cout <<"  "<<c<<": "<<mags[c-magmin]<<endl;
	cout <<endl;

	cout <<"Number of stars:      "<<numstars<<endl;
	cout <<"Number of galaxies:   "<<numgalaxies<<endl;
	cout <<"Brightest magnitude:  "<<magmin<<endl;
	cout <<"Dimmest magnitude:    "<<magmax<<endl;
	//cout <<"Color index range:   "<<colorindex_min<<"..."<<colorindex_max<<endl;



	 //copy data over to a gm image for output
	Image *stars=NULL;
	Image Stars;
	stars=&Stars;

	//stars->depth(8);
    //stars->magick("TIFF");
	//stars->matte(true);

	//char scratch[100];
	//sprintf(scratch,"%dx%d",width,height);
	//stars->size(scratch);
	//if (transparent) stars->read("xc:transparent");
	//else stars->read("xc:#00000000");

	stars->compressType(LZWCompression);
	stars->read(width,height,"RGBA",CharPixel,data);
    stars->magick("TIFF");


	cout <<"Writing to "<<filename<<"..."<<endl;
	stars->write(filename);



	cout << "Done!"<<endl;
	return 0;
}

/*! Return degrees from "246060". PGC uses this.
 */
double hms(const char *pos)
{
	if (*pos==' ') return 0;

	double v=(*pos-'0')*10;
	pos++;
	v+=(*pos-'0');
	pos++;

	double vv=(*pos-'0')*10;
	pos++;
	vv+=(*pos-'0');
	pos++;
	v+=vv/60.;

	vv=strtod(pos,NULL); //seconds has decimal
	v+=vv/3600.;

	v=v/24.*360;

	return v;
}

/*! Return degrees from "+906060". PGC uses this.
 */
double dms(const char *pos)
{
	if (*pos==' ') return 0;

	int sgn=1;
	if (*pos=='+') pos++;
	else if (*pos=='-') { sgn=-1; pos++; }

	 //degree portion, assume 0..90
	double v=(*pos-'0')*10;
	pos++;
	v+=(*pos-'0');
	pos++;

	double vv=(*pos-'0')*10;
	pos++;
	vv+=(*pos-'0');
	pos++;
	v+=vv/60.;

	vv=(*pos-'0')*10; //seconds has NO decimal here
	pos++;
	vv+=(*pos-'0');
	pos++;
	v+=vv/3600.;

	return v*sgn;
}

/*! Add color*alpha to pixel at x,y. If pixel has alpha, then make more opaque.
 */
void blendPixel(int x,int y, ColorRGB &color)
{
	if (x<0 || x>=width || y<0 || y>=height) return;

	int i=y*width*4+x*4;
	int r=data[i];
	int g=data[i+1];
	int b=data[i+2];
	int a=data[i+3];

	double aaa=color.alpha();
	if (alphaamp) {
		aaa+=alphaamp; if (aaa>1) aaa=1;
	}
	int aa=aaa*255;
	int rr=aaa*color.red()*255;
	int gg=aaa*color.green()*255;
	int bb=aaa*color.blue()*255;
	//int rr=color.red()*255;
	//int gg=color.green()*255;
	//int bb=color.blue()*255;

	r+=rr; if (r>255) r=255;
	g+=gg; if (g>255) g=255;
	b+=bb; if (b>255) b=255;
	a+=aa; if (a>255) a=255;

	data[i  ]=r;
	data[i+1]=g;
	data[i+2]=b;
	data[i+3]=a;
}

void point(double x,double y, ColorRGB &color, double span)
{
	if (span<=1) {
		 //near enough to equator that it is just a single pixel
		blendPixel(x,y,color);

	} else if (x-span/2<0) {
		 //draw partial, wraps around 
		for (int c=0; c<x+span/2; c++) blendPixel(c,y,color); //draw line 0 to x+span/2
		for (int c=x+width-span/2; c<width; c++) blendPixel(c,y,color); //draw line x+width-span/2, to width

	} else if (x+span/2>width) {
		 //wraps around to beginning
		for (int c=0; c<x-width+span/2; c++) blendPixel(c,y,color); //line 0 to x-width+span/2
		for (int c=x-span/2; c<width; c++) blendPixel(c,y,color); //line x-span/2 to width

	} else {
		 //fits without wraparound
		for (int c=x-span/2; c<x+span/2; c++) blendPixel(c,y,color); //line x-span/2 to width
	}
}

void dataEllipse(int xp,int yp, double xr,double yr, ColorRGB &color)
{
	double xspan;

	if (usehalo) {
		 //may halo reference to data
		 //rendered on in a very naive way, just rectangular copy, no additional span correction

		double oldamp=alphaamp;
		alphaamp=0;

		ColorRGB col;
		int i;
		double a,r,g,b;
		int sx,sy;
		int hx,hy;

		for (int y=yp-yr; y<yp+yr; y++) {
		  sy=y-(yp-yr);
		  for (int x=xp-xr; x<xp+xr; x++) {
		    sx=x-(xp-xr);

			r=color.red();
			g=color.green();
			b=color.blue();

			hy=sy/(2.*yr)*halowidth;
			hx=sx/(2.*xr)*halowidth;
			i=hy*halowidth + hx;
			a=halo[i]/255.;

			//cerr <<"hx,hy:"<<hx<<','<<hy<<"  i,a:"<<i<<"  "<<a<<endl;

			//col.alpha(a*color.alpha());
			col.alpha(a * (1-.25*color.alpha()));
			if (a>.5) {
				r=r*(1-a) + 1*(a);
				g=g*(1-a) + 1*(a);
				b=b*(1-a) + 1*(a);
			}
			col.red  (r);
			col.green(g);
			col.blue (b);

			blendPixel(x,y, col);
		  }
		  //cerr <<endl;
		}
		
		alphaamp=oldamp;

	} else {
		 //draw (alas non-antialiased) pixels in a circle
		for (int y=yp-yr; y<yp+yr; y++) {
		  xspan=xr*sqrt(1-(y-yp)*(y-yp)/yr/yr);
		  for (int x=xp-xspan; x<xp+xspan; x++) {
			blendPixel(x,y, color);
		  }
		}
	}
}


/*! span is how much horizontally to stretch out one pixel.
 * Vertical is not stretched. Note this is not accurate.. ellipses do not actually project into ellipses.
 */
void ellipse(double x,double y, double xs,double ys, ColorRGB &color, double span)
{
	 //tinker with colors to make more pleasing
	////color.alpha(0);
	//color.alpha(color.alpha()*.25);
	////color.red(color.red()*.5);

	//if (color.blue()<.4) color.blue(color.red()*.75);
	//if (color.green()<.4) color.green(color.red()*.75);

//	if (usehalo) {
//		xs*=usehalo;
//		ys*=usehalo;
//	}

	 //xs and ys were diameters, convert to radii
	xs*=span/2; //expand when necessary
	ys/=2;


	if (x-xs<0) {
		 //near 0, draw partial, wraps around to far
		dataEllipse(x,y, xs,ys, color);
		dataEllipse(x+width,y, xs,ys, color);

	} else if (x+xs>width) {
		 //near far edge, wraps around to begining
		dataEllipse(x-width,y, xs,ys, color);
		dataEllipse(x,y, xs,ys, color);

	} else {
		 //fits without wraparound
		dataEllipse(x,y, xs,ys, color);
	}
}







/**
 * This draws a pixel or small circle to the screen.
 */
void drawStar(float ra, float dec, float vmag, float bmag)
{
  float x, y;
  ColorRGB c;
  float index = bmag-vmag;
  indexToRgb(index, vmag, c);  // Get a color for the star
  //stroke(c);
  //fill(c);
  
  // map to screen coordinates
  x = map(ra, 0, 360, 0, width);
  y = map(dec, -90, 90, height, 0);

  //cout <<"x:"<<x<<"  y:"<<y<<endl;
  
  
  // if Galactic Coordinates..
  if (galactic) {
  	  //cout <<"ra:"<<radians(ra)<<"  dec:"<<radians(dec)<<endl;

	  flatvector galacticCoord = Eq2Gal(radians(ra), radians(dec));
	  x = galacticCoord.x;
	  
	  // Put the center of the Galaxy at the center of the image
	  if(x > 180)
		x -= 360;
	  
	  // map to screen coordinates
	  x = map(x, -180, 180, 0, width);
	  y = map(galacticCoord.y, -90, 90, height, 0);
  }
  
   //figure out how much the star is stretched out in final image
  double span=1;
  if (galactic) dec=180.*(y-height/2)/height;
  double declination_radians=dec/180.*M_PI;
  if (declination_radians==M_PI) span=width;
  else span=1/cos(declination_radians);
  if (span>width) span=width;


  // For bright stars draw a cicle, otherwise just a pixel
  if (vmag < bigthreshhold)
  {
    float s = map(vmag, bigthreshhold, -1, 2, maxstarsize);
	ellipse(x, y, s, s,  c, span);
  }
  else
  {
	point(x, y, c, span);
  }
}


//! Figure out how much one normal pixel is stretched out horizontally torward the poles.
double GetSpan(int y)
{
  double span=1;
  double dec=180.*(y-height/2)/height;
  double declination_radians=dec/180.*M_PI;
  if (declination_radians==M_PI) span=width;
  else span=1/cos(declination_radians);
  if (span>width) span=width;
  return span;
}


//for use in Eq2Gal:
float an          = radians(32.93192);   // Galactic long of asc node on equator
float ngpRa       = radians(192.85948);  // RA of North Galactic Pole
float ngpDec      = radians(27.12825);   // Dec of North Galactic Pole
float cos_ngpDec  = cos(ngpDec);
float sin_ngpDec  = sin(ngpDec);
static float SMALL = 1e-20;

/**
 * Convert Equtorial Coordinates to Galactic Coordinates
 * Based on code from libastro. You are not expected to understand this.
 */
flatvector Eq2Gal(float ra, float dec)
{
  float sin_dec, cos_dec, a, cos_a, sin_a, b, square, c, d; 
  float lat_gal;
  float lon_gal;
  
  cos_dec = cos(dec);
  sin_dec = sin(dec);
  a = ra - ngpRa;
  cos_a = cos(a);
  sin_a = sin(a);
  b = cos_a;
  square = (cos_dec*cos_ngpDec*b) + (sin_dec*sin_ngpDec);

  // Galactic Latitude
  lat_gal = asin(square);

  c = sin_dec - (square*sin_ngpDec);
  d = cos_dec*sin_a*cos_ngpDec;
  if (abs(d) < SMALL)
    d = SMALL;

  // Galactic Longitude
  lon_gal = atan(c/d) + an;
    
  if (d < 0) lon_gal += PI;
  if (lon_gal < 0) lon_gal += TWO_PI;
  if (lon_gal > TWO_PI) lon_gal -= TWO_PI;
  
  return flatvector(degrees(lon_gal), degrees(lat_gal));
}





/**
 * Method for mapping color index to a Color object
 * Based on http://www.vendian.org/mncharity/dir3/starcolor/
 * numbers are more of less made up but match
 * vendian.org's results pretty well.
 *
 * Just in case the website fades, mncharity's results:
 * 
 *   B-V      Teff              B-V     Teff              B-V     Teff             B-V      Teff  
 * -0.40    113017   #9bb2ff   0.25     7483   #eeefff   0.90     5052   #ffe8ce   1.55     3892   #ffd29c
 * -0.35     56701   #9eb5ff   0.30     7218   #f3f2ff   0.95     4948   #ffe6ca   1.60     3779   #ffd096
 * -0.30     33605   #a3b9ff   0.35     6967   #f8f6ff   1.00     4849   #ffe5c6   1.65     3640   #ffcc8f
 * -0.25     22695   #aabfff   0.40     6728   #fef9ff   1.05     4755   #ffe3c3   1.70     3463   #ffc885
 * -0.20     16954   #b2c5ff   0.45     6500   #fff9fb   1.10     4664   #ffe2bf   1.75     3234   #ffc178
 * -0.15     13674   #bbccff   0.50     6285   #fff7f5   1.15     4576   #ffe0bb   1.80     2942   #ffb765
 * -0.10     11677   #c4d2ff   0.55     6082   #fff5ef   1.20     4489   #ffdfb8   1.85     2579   #ffa94b
 * -0.05     10395   #ccd8ff   0.60     5895   #fff3ea   1.25     4405   #ffddb4   1.90     2150   #ff9523
 * -0.00      9531   #d3ddff   0.65     5722   #fff1e5   1.30     4322   #ffdbb0   1.95     1675   #ff7b00
 *  0.05      8917   #dae2ff   0.70     5563   #ffefe0   1.35     4241   #ffdaad   2.00     1195   #ff5200
 *  0.10      8455   #dfe5ff   0.75     5418   #ffeddb   1.40     4159   #ffd8a9
 *  0.15      8084   #e4e9ff   0.80     5286   #ffebd6   1.45     4076   #ffd6a5
 *  0.20      7767   #e9ecff   0.85     5164   #ffe9d2   1.50     3989   #ffd5a1 
 *
 * The vmag value is mapped to the alpha of color.
 */
void indexToRgb(float index, float vmag, ColorRGB &color)
{
  //float bright = map(vmag,  magnitude, -1,  0, 1); //map magnitude range to 0..1
  float bright = map(vmag,  magnitude, bigthreshhold+.0001,  0, 1); //map magnitude range of point starsto 0..1
 
  //if (vmag>bigthreshhold) cout <<"bright:"<<bright<<endl;

  double r = indexToRed(index);
  double g = indexToGreen(index);
  double b = indexToBlue(index);
  
  color.red(r);
  color.green(g);
  color.blue(b);
  color.alpha(bright);
  //cerr <<"bright: "<<bright<<endl;
}


double indexToRed(float index)
{  
  if (index<.5) {
	double r=110*index + 200;
	r=r/255;
	r=constrain(r, 0,1);
    return r;

  } else return 1;
}

double indexToGreen(float index)
{
  if (index < 0.4)
  {
    double g = (94.4*index + 217.2)/255;
    g = constrain(g, 0, 1);
    return g;

  }
  else if (index <= 1.7)
  {
    double g = (-42.3*index + 272)/255;
    g = constrain(g, 0, 1);
    return g;

  }
  else
  {
    double g = (-667*index + 1333)/255;
    g = constrain(g, 0, 1);
    return g;

  }
}

double indexToBlue(float index)
{  
  if (index < .5) {
	return 1;

  } else if (index<1.7) {
	double bl=-200*index + 300;
	bl=bl/255;
	bl=constrain(bl, 0,1);
	return bl;

  } else {
    float bl = (-500*index + 1000)/255;

    bl= constrain(bl, 0, 1);
    return bl;
  }
}


/*! Create a square w x h grayscale image (no alpha), where 0 is transparent, 1 is full opaque.
 * halo size is the ratio of halo diameter to star diameter.
 * These are copied to star map for bigger stars.
 *
 * w must be an even number.
 *
 * For debugging, writes out to halo.png.
 */
void CreateStockHalo(int w,double halosize, unsigned char *halodata)
{
	double r=w/2/halosize; //radius of main star
	double r2=r*r;         //square of radius of main star
	int w2=w*w/4;          //square of halo radius
	int cx=w/2, cy=w/2;    //center of halo
	double d2;             //temp, distance of current point to center
	double hstart=1;       //what radius to start full opaque halo
	hstart*=r;

	double v;

	for (int y=0; y<w; y++) {
	  for (int x=0; x<w; x++) {
		d2=(x-cx)*(x-cx)+(y-cy)*(y-cy); //distance squared of current pixel to center

		if (d2<r2) { //within main star
			halodata[y*w+x]=255;//totally opaque

		} else if (d2<w2) { //within halo
			//draw halo pixel
			//v: 1==opaque, 0==transparent
			v=((w/2) - sqrt(d2)) / (w/2-hstart);

			//halo dropoff..
			//1. so far it is linear

			//2. maybe v^2
			//v=v*v;// *** testing different halo types

			//3. sine wave
			v=constrain(v, 0,1);
			v=.5+.5*sin((v-.5)*M_PI);

			halodata[y*w+x]=constrain(v*255, 0,255);

		} else  {
			//totally transparent, outside of halo
			halodata[y*w+x]=0;
		}

	  }
	}


	Image haloimage;
	haloimage.read(halowidth,halowidth,"I",CharPixel,halodata);
    haloimage.magick("PNG");
	haloimage.write("halo.png");
}



