


#include "catalogs.h"

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include </usr/include/GraphicsMagick/Magick++.h>

#include <lax/strmanip.h>

#include <lax/lists.cc>


using namespace std;
using namespace Laxkit;
using namespace Magick;


#define PI (M_PI)
#define TWO_PI (2*M_PI)

//map: scale from [omin,omax] to [nmin,nmax]
#define map(value, omin,omax, nmin,nmax) (nmin+(value-omin)/(omax-omin)*(nmax-nmin))
#define constrain(v, min,max) (v<min ? min : (v>max ? max : v))
#define radians(deg) (deg/180.*M_PI)
#define degrees(rad) (rad*180./M_PI)

 //integer color average
#define AVG(a,b, r) ((a*r+b*(255-r))/255)


//------------------------------- RenderContext -----------------------------------
/*! \class RenderContext
 * Info about how, what, and where to render.
 */



RenderContext::RenderContext()
{
	 //curve objects:
	ramp=NULL;
	blowout=NULL;
	index_r=NULL;
	index_g=NULL;
	index_b=NULL;


	projectfile=NULL;
	filename   =newstr("stars.tif");
	width      =8192;
	height     =4096;
	galactic   =0;
	transparent=0;

	maximum_magnitude=-5;
	minimum_magnitude=15.1;

	maxstarsize  =3; //pixels of biggest star
	bigthreshhold=4; //magnitude < this get big treatment
	alphaamp     =0;

	usehalo  =5;
	halo     =NULL;
	halowidth=0;

	data=NULL;
	catalog=NULL;
}


//------------------------------- Random Catalog for previewing mainly -----------------------------------


/*! \class RandomCatalog
 * Collection of stars in memory for previewing.
 */
RandomCatalog::RandomCatalog(const char *nname, int num)
  : Catalog(nname,NULL,RandomMemory)
{
	if (num<=0) num=1000;
	numpoints=num;
	color_index=new double[numpoints];
	color_mag  =new double[numpoints];
	asc        =new double[numpoints];
	dec        =new double[numpoints];

	Repopulate(numpoints,0);
}

//! Convert cartesian space to spherical surface. (radians)
void rect_to_sphere(double x,double y,double z, double *asc,double *dec)
{
	*asc=atan2(y,x);
	*dec=atan(z/sqrt(x*x+y*y));
}

int RandomCatalog::Repopulate(int num, int spherical)
{
	if (num>numpoints) {
		numpoints=num;
		delete[] color_index;
		delete[] color_mag;
		delete[] asc;
		delete[] dec;
		color_index=new double[numpoints];
		color_mag  =new double[numpoints];
		asc        =new double[numpoints];
		dec        =new double[numpoints];
	}

	double x,y,z;
	for (int c=0; c<numpoints; c++) {
		color_index[c]=drand48()*(2.5)-.5;
		color_mag[c]  =drand48()*(15);
		
		 //create random point in unit cube, map to sphere..
		 //creates even distribution compared to just random asc/dec
		if (spherical) {
			x=drand48();
			y=drand48();
			z=drand48();

			if (x==0 && y==0 && z==0) x=1;
			rect_to_sphere(x,y,z, &asc[c], &dec[c]);
			asc[c]/=2*M_PI; //convert to unit range
			if (asc[c]<0) asc[c]+=.5;
			dec[c]/=M_PI/2; //convert to unit range
			if (dec[c]<0) dec[c]+=.5;

		} else {
			 //for simple rectilinear:
			asc[c]=drand48();
			dec[c]=drand48();
		}
	}

	return numpoints;
}

int RandomCatalog::Render(RenderContext *context)
{
	return 1;
}

int RandomCatalog::Render(RenderContext *context, unsigned char *data,int ww,int hh)
{
    double Ra;   //right ascension (longitude)
    double Dec;  //declination (latitude)
    double vmag; //visual mag
    double index;//blue mag


	 //preempt to draw on provided surface
	unsigned char *olddata=context->data;
	int oldw=context->width;
	int oldh=context->height;

	context->data  =data;
	context->width =ww;
	context->height=hh;

	for (int c=0; c<numpoints; c++) {

		vmag  = color_mag[c];
		index = color_index[c];
		Ra    = asc[c];
	 	Dec   = dec[c];

		if (vmag>context->minimum_magnitude || vmag<context->maximum_magnitude) continue;

		drawStarSimple(context, Ra, Dec, index, vmag);
	}


	 //restore actual surface
	context->data  =olddata;
	context->width =oldw;
	context->height=oldh;  

	return 0;
}


//------------------------------- Catalog -----------------------------------

/*! \class Catalog
 * Hold information about a sky catalog,
 * such as the Tycho 2 Star Catalog, or the Principal Galaxy Catalog.
 */

Catalog::Catalog(const char *nname, const char *nfile, CatalogTypes ntype)
{
	name=newstr(nname);
	filename=newstr(nfile);
	type=ntype;

	min_mag_cutoff=max_mag_cutoff=0;
	
	 //stats stuff:
	magnitude_distribution=NULL;
	nummags=0;
	minimum_magnitude=maximum_magnitude=0;
}

Catalog::~Catalog()
{
	if (name) delete[] name;
	if (filename) delete[] filename;
}


//! Default is to just fopen(filename), return nonzero for error, or 0 for success.
int Catalog::OpenCatalog()
{
	return 1;

//	FILE *f=NULL;
//	if (f) {
//		cerr <<"Catalog already open! Closing before proceeding..."<<endl;
//		CloseCatalog();
//	}
//
//	f=fopen(filename,"r");
//	if (!f) {
//		cerr <<" --Fail!-- Could not open "<<(name?name:"catalog")<<": "<<filename<<endl;
//		return 1;
//	}
//	return 0;
}

//! Default is to just fclose().
int Catalog::CloseCatalog()
{
	//if (!f) return 1;
	//fclose(f);
	return 0;
}




int Catalog::Render(RenderContext *context)
{
	if (type==PrincipalGalaxy) return Process_PGC(context);
	if (type==Tycho2)          return Process_Tycho(context);

	return 1;
}




//------------------------- Principal Galaxy Catalog processing ---------------------------
int Process_PGC(RenderContext *rr)
{
	const char *pgc_file=rr->catalog->filename;


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


	int numgalaxies=0;
	if (pgc_file) {
		char *line=NULL;
		size_t n=0;
		ssize_t c;
		FILE *f=fopen(pgc_file,"r");
		if (!f) {
			cerr <<" --Fail!-- Could not open Principal Galaxy Catalog: "<<pgc_file<<endl;
			return 1;
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
		 
			if (vmag>rr->minimum_magnitude || vmag<rr->maximum_magnitude) continue;
			if (Ra==0 && Dec==0) continue;

			if (vmag<magmin) magmin=vmag;
			if (vmag>magmax) magmax=vmag;
			mags[int(vmag-magmin)]++;


			drawStar(rr, Ra, Dec, vmag, bmag);

			numgalaxies++;
			if (numgalaxies%10000==0) cout <<"+\n";

		} while (!feof(f));
		if (line) free(line);

		fclose(f);
	}

	return 0;
}



//-------------------------- Tycho star catalog processing ---------------------------------

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

int Process_Tycho(RenderContext *rr)
{
	const char *tycho_file=rr->catalog->filename;



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



	if (tycho_file) {
		char *line=NULL;
		size_t n=0;
		ssize_t c;
		FILE *f=fopen(tycho_file,"r");
		if (!f) {
			cerr <<" --Fail!-- Could not open Tycho catalog: "<<tycho_file<<endl;
			return (1);
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
		 
			if (vmag>rr->minimum_magnitude || vmag<rr->maximum_magnitude) continue;
			if (Ra==0 && Dec==0) continue;

			if (vmag<magmin) magmin=vmag;
			if (vmag>magmax) magmax=vmag;
			mags[int(vmag-magmin)]++;

			 //find stats about colors:
			//colorindex=bmag-vmag;
			//if (colorindex<colorindex_min) colorindex_min=colorindex;
			//if (colorindex>colorindex_max) colorindex_max=colorindex;


			drawStar(rr, Ra, Dec, vmag, bmag);

			numstars++;
			if (numstars%100000==0) cout <<".\n";

		} while (!feof(f));
		if (line) free(line);

		fclose(f);
	}

	return 0;
}



//------------------------------------ misc file parsing functions ---------------------------

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




//------------------------------------ Rendering code functions ---------------------------


/*! Add color*alpha to pixel at x,y. If pixel has alpha, then make more opaque.
 */
void blendPixel(RenderContext *context, int x,int y, StarColor &color)
{
	if (x<0 || x>=context->width || y<0 || y>=context->height) return;

	int i=y*context->width*4+x*4;
	int r=context->data[i];
	int g=context->data[i+1];
	int b=context->data[i+2];
	int a=context->data[i+3];

	double aaa=color.alpha();
	if (context->alphaamp) {
		aaa+=context->alphaamp; if (aaa>1) aaa=1;
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

	//context->data[i  ]=r;
	//context->data[i+1]=g;
	//context->data[i+2]=b;
	//context->data[i+3]=a;

	context->data[i  ]=r;
	context->data[i+1]=g;
	context->data[i+2]=b;
	context->data[i+3]=255;

	cerr <<"blend "<<x<<','<<y<<endl;
}

void point(RenderContext *rr, double x,double y, StarColor &color, double span)
{
	if (span<=1) {
		 //near enough to equator that it is just a single pixel
		blendPixel(rr, x,y,color);

	} else if (x-span/2<0) {
		 //draw partial, wraps around 
		for (int c=0; c<x+span/2; c++) blendPixel(rr, c,y,color); //draw line 0 to x+span/2
		for (int c=x+rr->width-span/2; c<rr->width; c++) blendPixel(rr, c,y,color); //draw line x+width-span/2, to width

	} else if (x+span/2>rr->width) {
		 //wraps around to beginning
		for (int c=0; c<x-rr->width+span/2; c++) blendPixel(rr, c,y,color); //line 0 to x-width+span/2
		for (int c=x-span/2; c<rr->width; c++) blendPixel(rr, c,y,color); //line x-span/2 to width

	} else {
		 //fits without wraparound
		for (int c=x-span/2; c<x+span/2; c++) blendPixel(rr, c,y,color); //line x-span/2 to width
	}
}

void dataEllipse(RenderContext *rr, int xp,int yp, double xr,double yr, StarColor &color)
{
	double xspan;

	if (rr->usehalo) {
		 //may halo reference to data
		 //rendered on in a very naive way, just rectangular copy, no additional span correction

		double oldamp=rr->alphaamp;
		rr->alphaamp=0;

		StarColor col;
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

			hy=sy/(2.*yr)*rr->halowidth;
			hx=sx/(2.*xr)*rr->halowidth;
			i=hy*rr->halowidth + hx;
			a=rr->halo[i]/255.;

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

			blendPixel(rr, x,y, col);
		  }
		  //cerr <<endl;
		}
		
		rr->alphaamp=oldamp;

	} else {
		 //draw (alas non-antialiased) pixels in a circle
		for (int y=yp-yr; y<yp+yr; y++) {
		  xspan=xr*sqrt(1-(y-yp)*(y-yp)/yr/yr);
		  for (int x=xp-xspan; x<xp+xspan; x++) {
			blendPixel(rr, x,y, color);
		  }
		}
	}
}


/*! span is how much horizontally to stretch out one pixel.
 * Vertical is not stretched. Note this is not accurate.. ellipses do not actually project into ellipses.
 */
void ellipse(RenderContext *rr, double x,double y, double xs,double ys, StarColor &color, double span)
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
		dataEllipse(rr, x,y, xs,ys, color);
		dataEllipse(rr, x+rr->width,y, xs,ys, color);

	} else if (x+xs>rr->width) {
		 //near far edge, wraps around to begining
		dataEllipse(rr, x-rr->width,y, xs,ys, color);
		dataEllipse(rr, x,y, xs,ys, color);

	} else {
		 //fits without wraparound
		dataEllipse(rr, x,y, xs,ys, color);
	}
}





/**
 * This draws a pixel or small circle to the screen, in a simple rectilinear way.
 * ra and dec are assumed to be in range 0..1.
 */
void drawStarSimple(RenderContext *rr, double ra, double dec, double index, double vmag)
{
  double x, y;
  StarColor c;
  indexToRgb(rr, index, vmag, c);  // Get a color for the star
  //stroke(c);
  //fill(c);
  
  // map to screen coordinates
  x = map(ra,  0, 1,     0,   rr->width);
  y = map(dec, 0, 1, rr->height, 0);

  //cout <<"x:"<<x<<"  y:"<<y<<endl;
  
  
   //figure out how much the star is stretched out in final image
  double span=1;


  // For bright stars draw a cicle, otherwise just a pixel
  if (vmag < rr->bigthreshhold)
  {
    float s = map(vmag, rr->bigthreshhold, -1, 2, rr->maxstarsize);
	ellipse(rr, x, y, s, s,  c, span);
  }
  else
  {
	point(rr, x, y, c, span);
  }
}



/**
 * This draws a pixel or small circle to the screen.
 */
void drawStar(RenderContext *rr, float ra, float dec, float vmag, float bmag)
{
  float x, y;
  StarColor c;
  float index = bmag-vmag;
  indexToRgb(rr, index, vmag, c);  // Get a color for the star
  //stroke(c);
  //fill(c);
  
  // map to screen coordinates
  x = map(ra, 0, 360, 0, rr->width);
  y = map(dec, -90, 90, rr->height, 0);

  //cout <<"x:"<<x<<"  y:"<<y<<endl;
  
  
  // if Galactic Coordinates..
  if (rr->galactic) {
  	  //cout <<"ra:"<<radians(ra)<<"  dec:"<<radians(dec)<<endl;

	  flatvector galacticCoord = Eq2Gal(radians(ra), radians(dec));
	  x = galacticCoord.x;
	  
	  // Put the center of the Galaxy at the center of the image
	  if(x > 180)
		x -= 360;
	  
	  // map to screen coordinates
	  x = map(x, -180, 180, 0, rr->width);
	  y = map(galacticCoord.y, -90, 90, rr->height, 0);
  }
  
   //figure out how much the star is stretched out in final image
  double span=1;
  if (rr->galactic) dec=180.*(y-rr->height/2)/rr->height;
  double declination_radians=dec/180.*M_PI;
  if (declination_radians==M_PI) span=rr->width;
  else span=1/cos(declination_radians);
  if (span>rr->width) span=rr->width;


  // For bright stars draw a cicle, otherwise just a pixel
  if (vmag < rr->bigthreshhold)
  {
    float s = map(vmag, rr->bigthreshhold, -1, 2, rr->maxstarsize);
	ellipse(rr, x, y, s, s,  c, span);
  }
  else
  {
	point(rr, x, y, c, span);
  }
}


//! Figure out how much one normal pixel is stretched out horizontally torward the poles.
double GetSpan(int y, int width,int height)
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



/*! Create a square w x h grayscale image (no alpha), where 0 is transparent, 1 is full opaque.
 * halo size is the ratio of halo diameter to star diameter.
 * These are copied to star map for bigger stars.
 *
 * w must be an even number.
 *
 * For debugging, writes out to halo.png.
 */
void CreateStockHalo(int w,double halosize, unsigned char *halodata,
					 const char *format, Laxkit::CurveInfo *ramp, Laxkit::CurveInfo *blowout)
{
	double r=w/2/halosize; //radius of main star
	double r2=r*r;         //square of radius of main star
	int w2=w*w/4;          //square of halo radius
	int cx=w/2, cy=w/2;    //center of halo
	double d2;             //temp, distance of current point to center
	double hstart=1;       //what radius to start full opaque halo
	hstart*=r;

	int i;
	double v;
	int vv;
	int usecolor=(*format!='g');
	int col_b=255;
	int col_g=0;
	int col_r=0;

	for (int y=0; y<w; y++) {
	  for (int x=0; x<w; x++) {
		if (*format=='g') i=y*w+x; else i=4*(y*w+x);

		d2=(x-cx)*(x-cx)+(y-cy)*(y-cy); //distance squared of current pixel to center

		if (d2<r2) { //within main star
			if (usecolor) {
				vv=ramp->lookup[255];
				halodata[i+0]=AVG(255,col_b, blowout->lookup[vv]); //b
				halodata[i+1]=AVG(255,col_g, blowout->lookup[vv]); //g
				halodata[i+2]=AVG(255,col_r, blowout->lookup[vv]);//r
				halodata[i+3]=vv; //a
			} else {
				halodata[i]=ramp->lookup[255];//totally opaque
			}

		} else if (d2<w2) { //within halo
			//draw halo pixel
			//v: 1==opaque, 0==transparent
			v=((w/2) - sqrt(d2)) / (w/2-hstart);

			//map to sine wave ---> using curve ramp now
			//v=constrain(v, 0,1);
			//v=.5+.5*sin((v-.5)*M_PI);

			if (usecolor) {
				vv=ramp->lookup[constrain(int(v*255), 0,255)];
				halodata[i+0]=AVG(255,col_b, blowout->lookup[vv]); //b
				halodata[i+1]=AVG(255,col_g, blowout->lookup[vv]); //g
				halodata[i+2]=AVG(255,col_r, blowout->lookup[vv]);//r
				halodata[i+3]=vv;
				//halodata[i+3]=ramp->lookup[vv];
			} else {
				halodata[i]=ramp->lookup[constrain(int(v*255), 0,255)];
			}

		} else  {
			//totally transparent, outside of halo
			if (usecolor) {
				halodata[i+0]=0;
				halodata[i+1]=halodata[i+2]=halodata[i+3]=0;
			} else {
				halodata[i]=0;
			}
		}

	  }
	}


	//Image haloimage;
	//if (usecolor) haloimage.read(w,w,"ARGB",CharPixel,halodata);
	//else haloimage.read(w,w,"I",CharPixel,halodata);
    //haloimage.magick("PNG");
	//haloimage.write("halo.png");
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
 * The vmag value is mapped to the alpha of color for point stars.
 */
void indexToRgb(RenderContext *rr, float index, float vmag, StarColor &color)
{
  //float bright = map(vmag,  magnitude, -1,  0, 1); //map magnitude range to 0..1
  float bright = map(vmag,  rr->minimum_magnitude, rr->bigthreshhold+.0001,  0, 1); //map magnitude range of point starsto 0..1
 
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



//----------------------------- Actual Render and Export Function --------------------


int Render(RenderContext *context)
{

	 //create halo image if necessary
	if (!context->halo) {
		context->halowidth=context->maxstarsize*context->usehalo;
		context->halo=new unsigned char[context->halowidth*context->halowidth];
		CreateStockHalo(context->halowidth, context->usehalo, context->halo, "g",
						context->ramp,context->blowout);
	}


	 //-----allocate star image data
	unsigned char *ddata=new unsigned char[(long)context->width*context->height*4];

	int stride=context->width*4;
	context->data=ddata;
	memset(context->data,0,context->width*context->height*4);
	if (!context->transparent) {
		for (int y=0; y<context->height; y++) {
		  for (int x=0; x<context->width; x++) {
			//data[y*stride+x*4+0]=255;
			//data[y*stride+x*4+1]=255;
			//data[y*stride+x*4+2]=255;
			context->data[y*stride+x*4+3]=255;//opaque
		  }
		}
	}


	 //Render the stack of catalogs
	for (int c=0; c<context->catalogs.n; c++) {
		context->catalog=context->catalogs.e[c];
		context->catalog->Render(context);
	}




	int mags[50];
	int magmin=1000, magmax=-1000;
	memset(mags,0,50*sizeof(int));
	int numstars=0;
	int numgalaxies=0;


	//-------------- print summary ------------------
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
	stars->read(context->width,context->height,"RGBA",CharPixel,context->data);
    stars->magick("TIFF");


	cout <<"Writing to "<<context->filename<<"..."<<endl;
	stars->write(context->filename);


	cout << "Done!"<<endl;

	return 0;
}



