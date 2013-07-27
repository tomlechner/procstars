// Star Catalog Image Generator
// By Tom Lechner, tomlechner.com
// 2013
//
// This program creates star globes from the wonderful Tycho2 star catalog,
// and also the Principal Galaxy Catalog.
//
// You can run with or with the gui.
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
//
// Fyi, the ISS takes up about 30 arc seconds.
// The moon takes up about 30 arc minutes.
// at 8192 px wide, 23 pixels make 1 degree. 1/3 of a pixel is one arc minute, 1/200 of a pixel is 1 arc second.
// at 32768,        91 pixels make 1 degree, 1.5 pixels make one arc minute. 
//
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





#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <stdint.h>

#include <lax/laxoptions.h>
#include <lax/rowframe.h>
#include <lax/laxutils.h>

#include </usr/include/GraphicsMagick/Magick++.h>

#include <lax/laxoptions.h>
#include <lax/anxapp.h>
#include <lax/curvewindow.h>
#include <lax/strmanip.h>

#include "catalogs.h"


#define PI (M_PI)
#define TWO_PI (2*M_PI)

//map: scale from [omin,omax] to [nmin,nmax]
#define map(value, omin,omax, nmin,nmax) (nmin+(value-omin)/(omax-omin)*(nmax-nmin))
#define constrain(v, min,max) (v<min ? min : (v>max ? max : v))
#define radians(deg) (deg/180.*M_PI)
#define degrees(rad) (rad*180./M_PI)



using namespace std;
using namespace Magick;
using namespace Laxkit;
using namespace LaxFiles;


//options:


//stuff to make computations easier:
//double pixelwidth=1./width;





//--------------------------------- IndexWindow ---------------------------------
/*! \class IndexWindow
 * Control window for adjusting map of star color index to rgb color.
 */
class IndexWindow : public Laxkit::RowFrame
{
  public:
	int numsamples;
	int maxsamplevalue;
	int *lookup_r;
	int *lookup_g;
	int *lookup_b;

	CurveWindow *rr,*gg,*bb;
	IndexWindow();
	virtual ~IndexWindow();
	virtual int init();
	virtual void Refresh();
	virtual int Event(const EventData *e,const char *mes);
	virtual void Reset(int which=~0);
};

IndexWindow::IndexWindow()
  : RowFrame(NULL,"Index Window","Index Window",ANXWIN_ESCAPABLE|ROWFRAME_ROWS, 0,0,600,300,0, NULL,0,NULL)
{
	lookup_r=lookup_g=lookup_b=NULL;
	numsamples=0;
	maxsamplevalue=255;
}

IndexWindow::~IndexWindow()
{
	if (lookup_r) delete[] lookup_r;
	if (lookup_g) delete[] lookup_g;
	if (lookup_b) delete[] lookup_b;
}

int IndexWindow::init()
{
	anXWindow *last=NULL;

	CurveWindow *win=NULL;
	
	AddSpacer(75,25,50,50,  200,100,200,50, -1);


	last=win=rr=new CurveWindow(NULL,"Red","Red",0, 5,5,500,500,0, last,object_id,"red",
						 			 "Red", "B-V",-0.5,2,  "r",0,255);
	rr->AddPoint(.4,255);
	AddWin(win,1, 200,100,200,50,0,  200,100,200,50,0, -1);


	last=win=gg=new CurveWindow(NULL,"Green","Green",0, 5,5,500,500,0, last,object_id,"green",
							"Green", "B-V",-0.5,2,  "g",0,255);
	gg->MovePoint(0, -.5,170);
	gg->MovePoint(1, 2,0);
	gg->AddPoint(.4,255);
	gg->AddPoint(1.7,200);
	AddWin(win,1, 200,100,200,50,0,  200,100,200,50,0, -1);


	last=win=bb=new CurveWindow(NULL,"Blue","Blue",0, 5,5,500,500,0, last,object_id,"blue",
							"Blue", "B-V",-0.5,2, "b",0,255);
	bb->MovePoint(0, -.5,255);
	bb->MovePoint(1, 2,0);
	bb->AddPoint(.5,255);
	bb->AddPoint(1.7,150);
	AddWin(win,1, 200,100,200,50,0,  200,100,200,50,0, -1);


	last->CloseControlLoop();
	Sync(1);
	return 0;
}

/*! Reset to default colors.
 * which is &1 for red, &2 for green, &4 for blue
 */
void IndexWindow::Reset(int which)
{
	if (which&1) {
		rr->Reset();
		rr->AddPoint(.4,255);
	}

	if (which&2) {
		gg->Reset();
		gg->MovePoint(0, -.5,170);
		gg->MovePoint(1, 2,0);
		gg->AddPoint(.4,255);
		gg->AddPoint(1.7,200);
	}

	if (which&4) {
		bb->Reset();
		bb->MovePoint(0, -.5,255);
		bb->MovePoint(1, 2,0);
		bb->AddPoint(.5,255);
		bb->AddPoint(1.7,150);
	}

	needtodraw=1;
}

void IndexWindow::Refresh()
{
	if (!needtodraw) return;
	if (arrangedstate!=1) Sync(0);

	int pad=10;
	int w=wholelist.e[0]->w()-2*pad;
	int h=wholelist.e[0]->h()-2*pad;
	double pos, r,g,b;

	for (int c=0; c<h; c++) {
		pos=(double)c/(h-1)*(2.4)-.4;
		r=rr->f(pos)/255;
		g=gg->f(pos)/255;
		b=bb->f(pos)/255;

		foreground_color(r,g,b);
		draw_line(this, pad,pad+c, pad+w,pad+c);
	}

	needtodraw=0;
}

int IndexWindow::Event(const EventData *e,const char *mes)
{
	const SimpleMessage *m=dynamic_cast<const SimpleMessage*>(e);
	if (!m) return anXWindow::Event(e,mes);

	needtodraw=1;
	return 0;
}




//--------------------------------- main() --------------------------------

LaxOptions options;

void InitOptions()
{
	options.HelpHeader( "Star Catalog Image Generator\n"
						"By Tom Lechner, 2013\n"
						"Version 0.000000001");
    options.UsageLine("procstars [options]");

	options.Add("width",        'w', 1, "Width of equirectangular texture image",                      0, "8192");
	options.Add("output-file",  'o', 1, "Filename for generated texture image (will save as tif)",     0, "file");
	options.Add("catalog-file", 'c', 1, "File containing the Tycho 2 Star Catalog",                    0, "catalog");
	options.Add("pgc-file",     'p', 1, "File containing the Principal Galaxy Catalog",                0, "catalog");
	options.Add("magnitude",    'm', 1, "Minimum brightness. stars in catalog are 0 (bright) "
			                            "to 15 (very dim). 6 is dimmest to human naked eye.",          0, "15.1");
	options.Add("maxstarsize",  's', 1, "Pixels wide of brightest star",0, "10");
	options.Add("maxstarangle", 'a', 1, "Arc minutes wide of brightest star",0, "20");
	options.Add("bigthreshhold",'S', 1, "Magnitude brighter than this gets special size treatment",    0, "4");
	options.Add("alpha-amp",    'A', 1, "Adjustment to artifically pump up dim stars. Between 0 and 1",0, "0");
	options.Add("halo",         'H', 1, "Render halos on bright stars, number is ratio of "
		                                "halo diameter to star diameter",                              0, "7");
	options.Add("transparent",  'T', 0, "Put stars on transparency rather than a black background",    0, NULL);
	options.Add("galactic",     'g', 0, "Make it so the Milky Way is horizontal",                      0, NULL);

	options.Add("help",         'h', 0, "Show help and exit",                                          0, NULL);
	options.Add("help-html",    'l', 0, "Output help in html form and exit",                           0, NULL);
	options.Add("version",      'v', 0, "Show version and exit",                                       0, NULL);
}


int main(int argc, char **argv)
{
	InitializeMagick(*argv);


	RenderContext rr;

	 //process command line options
	InitOptions();
	int cc,index;

	cc=options.Parse(argc,argv, &index);
    if (cc==-2) {
        cerr <<"Missing parameter for "<<argv[index]<<"!!"<<endl;
        exit(0);
    }
    if (cc==-1) {
        cerr <<"Unknown option "<<argv[index]<<"!!"<<endl;
        exit(0);
    }


	 //initialize option variables
	const char *tycho_file=NULL;
	const char *pgc_file=NULL;
	double ang=-1;


    LaxOption *o;
    for (o=options.start(); o; o=options.next()) {
        switch(o->chr()) {
            case 'l': {
				options.HelpHtml(stdout);
				exit(0);
			  }
            case 'h': {
				options.Help(stdout);
				exit(0);
			  }
            case 'v': {
				cout << options.HelpHeader()<<endl;
				exit(0);
			  }

            case 'w': {
				rr.width=strtol(o->arg(),NULL,10);
				if (rr.width<2) { cerr <<"Badth width value!"<<endl; exit(1); }
				rr.height=rr.width/2;
              } break;

            case 'o': {
				makestr(rr.filename,o->arg());
              } break;

            case 'p': {
				pgc_file=o->arg();
              } break;

            case 'c': {
				tycho_file=o->arg();
              } break;

            case 'm': {
				rr.minimum_magnitude=strtod(o->arg(),NULL);
              } break;

            case 'g': {
				rr.galactic=1;
              } break;

            case 'T': {
				rr.transparent=1;
              } break;

            case 'H': {
				rr.usehalo=strtod(o->arg(),NULL);
              } break;

            case 's': {
				rr.maxstarsize=strtod(o->arg(),NULL);
              } break;

            case 'a': {
				ang=strtod(o->arg(),NULL);
              } break;

            case 'S': {
				rr.bigthreshhold=strtod(o->arg(),NULL);
              } break;

            case 'A': {
				rr.alphaamp=strtod(o->arg(),NULL);
              } break;

		}
	}


	if (ang>0) {
		rr.maxstarsize=rr.width*ang/60./360.;
	}
	if (rr.usehalo>1) {
		rr.halowidth=rr.maxstarsize*rr.usehalo;
		rr.halo=new unsigned char[rr.halowidth*rr.halowidth];
		CreateStockHalo(rr.halowidth, rr.usehalo, rr.halo);
	}


	//pixelwidth=1./width;



	 //-------print options summary
	cout <<endl;
	cout << "      Tycho catalog: "<<(tycho_file?tycho_file:"none")<<endl;
	cout << "     Galaxy catalog: "<<(pgc_file?pgc_file:"none")<<endl;
	cout << "      Outputting to: "<<rr.filename<<endl;
	cout << "         Dimensions: "<<rr.width<<" x "<<rr.height<<endl;
	cout << "Brightness at least: "<<rr.minimum_magnitude<<endl;
	cout << "        Halos after: "<<rr.bigthreshhold<<endl;


	 //fill in up render context
	if (pgc_file)   rr.catalogs.push(new Catalog("Principal Galaxy Catalog", pgc_file,   PrincipalGalaxy),1);
	if (tycho_file) rr.catalogs.push(new Catalog("Tycho 2 Star Catalog",     tycho_file, Tycho2),         1);


	 //------set up windows...
	anXApp app;
	app.init(argc,argv);

	app.addwindow(new IndexWindow());
	

	 //...and off we go!
	app.run();


	cout <<"--- Shutting down... ----\n";
	app.close();
	
	cout <<"------ Bye! -------\n";
	return 0;

} //main()




void ExportStarImage(RenderContext *rr)
{




	 //-----allocate star image data
	unsigned char *ddata=new unsigned char[(long)rr->width*rr->height*4];

	int stride=rr->width*4;
	rr->data=ddata;
	memset(rr->data,0,rr->width*rr->height*4);
	if (!rr->transparent) {
		for (int y=0; y<rr->height; y++) {
		  for (int x=0; x<rr->width; x++) {
			//data[y*stride+x*4+0]=255;
			//data[y*stride+x*4+1]=255;
			//data[y*stride+x*4+2]=255;
			rr->data[y*stride+x*4+3]=255;//opaque
		  }
		}
	}





	//-------------- print summary ------------------


	int mags[50];
	int magmin=1000, magmax=-1000;
	memset(mags,0,50*sizeof(int));
	int numstars=0;
	int numgalaxies=0;


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
	stars->read(rr->width,rr->height,"RGBA",CharPixel,rr->data);
    stars->magick("TIFF");


	cout <<"Writing to "<<rr->filename<<"..."<<endl;
	stars->write(rr->filename);


	cout << "Done!"<<endl;
	return;
}









