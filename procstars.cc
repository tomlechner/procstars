// Tycho 2 Star Catalog Image Generator
// By Tom Lechner, tomlechner.com
// 2013
//
//
// The bits by Tom Lechner are released with the MIT license.
// There are other bits adapted from Nathan Bergey's code below,
// which he released CC Attribution 3.0.
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

//Usage:
//  ./process-tycho -S 6 -m 12 -s 12 -w 32768
//  ./process-tycho -S 6 -m 12 -s 6 -w 16384
//  ./process-tycho -c Tycho.cat -o stars.tif
//
//
//This program creates star globes from the wonderful Tycho2 star catalog.
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
//but after magnitude of about 12.1, you begin to see striations where the satellite did not scan.
//
//Todo:
//  blow out big stars, but color the halo
//  find a galaxy catalog: pgc
//  project billboard stuff like nebulae, LMC and SMC, M31, space stations, etc



#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
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
const char *tycho_file="Tycho.dat";
int width=8192;
int height=4096;
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



//forward decs
class flatvector
{
  public:
	double x,y;
	flatvector(double xx,double yy) { x=xx; y=yy; }
};
flatvector Eq2Gal(float ra, float dec);
void indexToRgb(float bmag, float vmag, ColorRGB &color);
double indexToRed(float index);
double indexToGreen(float index);
double indexToBlue(float index);
void drawStar(float ra, float dec, float vmag, float bmag);
void CreateStockHalo(int w,double halosize, unsigned char *halo);



Image *stars=NULL;


void help()
{
	cout << "Tycho 2 Star Catalog Image Generator"<<endl;
	cout << "By Tom Lechner"<<endl;
	cout << "Version 0.000000001"<<endl;
	cout << "\nprocess-tycho [options]\n\n";

	cout << "--width, -w 8192            Width of equirectangular texture image"<<endl;
	cout << "--output-file, -o file      Filename for generated texture image (will save as tif)"<<endl;
	cout << "--catalog-file, -c catalog  File containing the Tycho 2 Star Catalog"<<endl;
	cout << "--magnitude, -m 15.1        Minimum brightness. stars in catalog are 0 (bright)\n"
			"                              to 15 (very dim). 6 is dimmest to human naked eye."<<endl;
	cout << "--transparent, -t           Put stars on transparency rather than a black background"<<endl;
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
	Image Stars;
	stars=&Stars;

     // parse options
    static struct option long_options[] = {
            { "width",          1, 0, 'w' },
            { "output-file",    1, 0, 'o' },
            { "catalog-file",   1, 0, 'c' },
            { "magnitude",      1, 0, 'm' },
            { "maxstarsize",    1, 0, 's' },
            { "maxstarangle",   1, 0, 'a' },
            { "bigthreshhold",  1, 0, 'S' },
            { "galactic",       0, 0, 'g' },
            { "halo",           1, 0, 'H' },
            { "transparent",    0, 0, 't' },
            { "alpha-amp",      1, 0, 'A' },
			{ "version",        0, 0, 'v' },
            { "help",           0, 0, 'h' },
            { 0,0,0,0 }
        };
    int cc,index;

	double ang=-1;
    while (1) {
        cc=getopt_long(argc,argv,":w:o:c:m:s:S:a:A:H:gtvh",long_options,&index);
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

            case 'c': {
				tycho_file=optarg;
              } break;

            case 'm': {
				magnitude=strtod(optarg,NULL);
              } break;

            case 'g': {
				galactic=1;
              } break;

            case 't': {
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


	cout << "      Tycho catalog: "<<tycho_file<<endl;
	cout << "      Outputting to: "<<filename<<endl;
	cout << "         Dimensions: "<<width<<" x "<<height<<endl;
	cout << "Brightness at least: "<<magnitude<<endl;





	unsigned char *ddata=new unsigned char[width*height*4];

	int stride=width*4;
	data=ddata;
	 //init data
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


	//stars->depth(8);
    //stars->magick("TIFF");
	//stars->matte(true);

	//char scratch[100];
	//sprintf(scratch,"%dx%d",width,height);
	//stars->size(scratch);
	//if (transparent) stars->read("xc:transparent");
	//else stars->read("xc:#00000000");



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


	 //--------------- Tycho star catalog processing ---------------
	char *line=NULL;
	size_t n=0;
	ssize_t c;
	FILE *f=fopen(tycho_file,"r");
	if (!f) {
		cerr <<" --Fail!-- Could not open Tycho catalog: "<<tycho_file<<endl;
		exit(1);
	}
	int numstars=0;
	do {

		c=getline(&line,&n,f);
		if (c<0) break;
		if (c==0) continue;


		 //the catalog lines are delimited by '|' characters, but also arranged on strict byte widths
        Ra   = strtod(line+15,NULL); //pieces[2],  mRA deglees, bytes 16-27
        Dec  = strtod(line+28,NULL); //pieces[3],  mDE degrees, bytes 29-40
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


	//--------------- Principal Galaxy Catalog processing -------------------
	//
	//



	//-------------- All done! print summary------------------

	 //statistics
	cout <<"\nNumber of stars:   "<<numstars<<endl;
	cout <<"Brightest magnitude: "<<magmin<<endl;
	cout <<"Dimmest magnitude:   "<<magmax<<endl;
	//cout <<"Color index range:   "<<colorindex_min<<"..."<<colorindex_max<<endl;
	cout <<"Magnitudes:"<<endl;
	for (int c=magmin; c<=magmax; c++) cout <<"  "<<c<<": "<<mags[c-magmin]<<endl;



	stars->compressType(LZWCompression);
	stars->read(width,height,"RGBA",CharPixel,data);
    stars->magick("TIFF");


	cout <<"Writing to "<<filename<<"..."<<endl;
	stars->write(filename);


	cout << "Done!"<<endl;
	return 0;
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
		 //wraps around to begining
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
		 //rendered on in a very naive way, just rectangular copy

		double oldamp=alphaamp;
		alphaamp=0;

		ColorRGB col;
		int i;
		double a;
		int sx,sy;
		int hx,hy;

		for (int y=yp-yr; y<yp+yr; y++) {
		  sy=y-(yp-yr);
		  for (int x=xp-xr; x<xp+xr; x++) {
		    sx=x-(xp-xr);

			col.red(color.red());
			col.green(color.green());
			col.blue(color.blue());

			hy=sy/(2.*yr)*halowidth;
			hx=sx/(2.*xr)*halowidth;
			i=hy*halowidth + hx;
			a=halo[i]/255.;

			//cerr <<"hx,hy:"<<hx<<','<<hy<<"  i,a:"<<i<<"  "<<a<<endl;

			col.alpha(a*color.alpha());
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
 * Vertical is not stretched. Not this is not accurate.. ellipses do not actually project into ellipses.
 */
void ellipse(double x,double y, double xs,double ys, ColorRGB &color, double span)
{
	 //tinker with colors to make more pleasing
	////color.alpha(0);
	//color.alpha(color.alpha()*.25);
	////color.red(color.red()*.5);

	//if (color.blue()<.4) color.blue(color.red()*.75);
	//if (color.green()<.4) color.green(color.red()*.75);

	if (usehalo) {
		xs*=usehalo;
		ys*=usehalo;
	}

	 //xs and ys were diameters
	xs*=span/2; //expand when necessary
	ys/=2;
	span=xs;


	if (x-xs/2<0) {
		 //near 0, draw partial, wraps around to far
		dataEllipse(x,y, xs,ys, color);
		dataEllipse(x+width,y, xs,ys, color);

	} else if (x+xs/2>width) {
		 //near far edge, wraps around to begining
		dataEllipse(x-width,y, xs,ys, color);
		dataEllipse(x,y, xs,ys, color);

	} else {
		 //fits without wraparound
		dataEllipse(x,y, xs,ys, color);
	}
}

/*! Create a square w x h grayscale image (no alpha), where 0 is transparent, 1 is full opaque.
 * halo size is the ratio of halo diameter to star diameter.
 * These are copied to star map for bigger stars.
 *
 * w must be an even number.
 */
void CreateStockHalo(int w,double halosize, unsigned char *halodata)
{
	double r=w/2/halosize; //radius of main star
	double r2=r*r;         //square of radius of main star
	int w2=w*w/4;          //square of halo radius
	int cx=w/2, cy=w/2;    //center of halo
	double d2;             //temp, distance of current point to center
	double hstart=0;       //what radius to start full opaque halo
	double v;

	for (int y=0; y<w; y++) {
	  for (int x=0; x<w; x++) {
		d2=(x-cx)*(x-cx)+(y-cy)*(y-cy); //distance squared of current pixel to center

		if (d2<r2) { //within main star
			halodata[y*w+x]=255;//totally opaque

		} else if (d2>w2) {
			halodata[y*w+x]=0;//totally transparent, outside of halo

		} else {
			//draw halo pixel
			v=((w/2)-sqrt(d2)-hstart)/(w/2-hstart);
			//v=v*v;// *** testing different halo types
			halodata[y*w+x]=constrain(v*255, 0,255);
		}
	  }
	}


	Image haloimage;
	haloimage.read(halowidth,halowidth,"I",CharPixel,halodata);
    haloimage.magick("PNG");
	haloimage.write("halo.png");
}


/**
 * The following is adapted from:
 *
 * A short piece of code to read the tycho.dat file which is the combined data 
 * from all the other files
 * Reads each line in the file and prints and image
 *
 * (c) Nathan Bergey 01-22-2010
 * Availible under Creative Commons Attribution 3.0 license.
 */

float an          = radians(32.93192);   // Galactic long of asc node on equator
float ngpRa       = radians(192.85948);  // RA of North Galactic Pole
float ngpDec      = radians(27.12825);   // Dec of North Galactic Pole
float cos_ngpDec  = cos(ngpDec);
float sin_ngpDec  = sin(ngpDec);
static float SMALL = 1e-20;
//PFont font;




/**
 * This draws a pixel or small circle to the screen.
 */
void drawStar(float ra, float dec, float vmag, float bmag)
{
  float x, y;
  ColorRGB c;
  indexToRgb(bmag, vmag, c);  // Get a color for the star
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

///**
// * Draws a grid and numbers (coordinates) on the screen
// */
//void printGrid()
//{
//  float x, y;
//  stroke(100);
//  strokeWeight(2);
//  fill(255);
//  
//  // x for Equtorial
//  for(int i = 0; i <= 24; i+= 2)
//  { 
//    x = map(i, 0, 24, 0, width);
//    line(x, 0, x, height);
//    text(i, x - 60, 50);
//  }
//
//  /* Uncomment for galactic coordinates
//   *
//  // x for Galactic
//  for(int i = -180; i <= 180; i+= 30)
//  { 
//    x = map(i, -180, 180, 0, width);
//    line(x, 0, x, height);
//    text(i, x - 60, 50);
//  }
//  */
// 
//  // y
//  for(int i = -90; i < 90; i += 30)
//  {
//    y = map(i, -90, 90, height, 0);
//    line(0, y, width, y);
//    text(i, 10, y + 50);
//  }
//}




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
 */
void indexToRgb(float bmag, float vmag, ColorRGB &color)
{
  float index = bmag - vmag;
  //float bright = map((vmag<bmag?vmag:bmag), 18, -1, 0, 1);
  float bright = map(vmag, 18, -1, 0, 1);
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

//--------------------------------older--------------------------------
double indexToRed_OLD(float index)
{  
  float a = 132.206;
  float b = 206.618;
  
  float r = ((a*index) + b)/255;
  r = constrain(r, 0, 1);
  
  return r;
}

double indexToGreen_OLD(float index)
{  
  if (index < 0.4)
  {
    float a = 92.9412;
    float b = 216.647;
    float g = ((a*index) + b)/255;
    g = constrain(g, 0, 1);
    return g;
  }
  else if (index <= 1.5)
  {
    float a = -33.7549;
    float b = 262.957;
    float g = ((a*index) + b)/255;
    g = constrain(g, 0, 1);
    return g;
  }
  else
  {
    float a = -564.678;
    float b = 210.636905;
    float c = 1.55;
    float g = (a*((index-c)*(index-c)) + b)/255;
    g = constrain(g, 0, 1);
    return g;
  }
}

double indexToBlue_OLD(float index)
{  
  if (index < 1.65)
  {
    float a = -84.1162;
    float b = 284.527;
    float bl = ((a*index) + b)/255;
    bl = constrain(bl, 0, 1);
    return bl;
  }
  else
  {
    float a = -530.0;
    float b = 1020.0;
    float bl = ((a*index) + b)/255;
    bl= constrain(bl, 0, 1);
    return bl;
  }
}

//-------------toms old:

/*! GraphicsMagick does not seem to apply alpha to pixels, so moving it down here juts for reference...
 */
void point_gm(double x,double y, ColorRGB &color, double span)
{
	if (span<=1) {
		 //near enough to equator that it is just a single pixel
		stars->pixelColor(x,y,color);

	} else if (x-span/2<0) {
		 //draw partial, wraps around 
		for (int c=0; c<x+span/2; c++) stars->pixelColor(c,y,color); //draw line 0 to x+span/2
		for (int c=x+width-span/2; c<width; c++) stars->pixelColor(c,y,color); //draw line x+width-span/2, to width

	} else if (x+span/2>width) {
		 //wraps around to begining
		for (int c=0; c<x-width+span/2; c++) stars->pixelColor(c,y,color); //line 0 to x-width+span/2
		for (int c=x-span/2; c<width; c++) stars->pixelColor(c,y,color); //line x-span/2 to width

	} else {
		 //fits without wraparound
		for (int c=x-span/2; c<x+span/2; c++) stars->pixelColor(c,y,color); //line x-span/2 to width
	}
}

void ellipse_gm(double x,double y, double xs,double ys, ColorRGB &color, double span)
{
	 //tinker with colors to make more pleasing
	//color.alpha(0);
	color.alpha(color.alpha()*.25);
	//color.red(color.red()*.5);
	if (color.blue()<.4) color.blue(color.red()*.95);
	if (color.green()<.4) color.green(color.red()*.95);


	xs*=span; //expand when necessary
	stars->fillColor(color);
	stars->strokeWidth(0);
	stars->strokeColor(color);

	span=xs;


	if (x-xs/2<0) {
		 //near 0, draw partial, wraps around to far
		stars->draw(DrawableEllipse(x,y, xs,ys, 0,360));
		stars->draw(DrawableEllipse(x+width,y, xs,ys, 0,360));

	} else if (x+xs/2>width) {
		 //near far edge, wraps around to begining
		stars->draw(DrawableEllipse(x-width,y, xs,ys, 0,360));
		stars->draw(DrawableEllipse(x,y, xs,ys, 0,360));

	} else {
		 //fits without wraparound
		stars->draw(DrawableEllipse(x,y, xs,ys, 0,360));
	}
}

