//  Procstars, Star catalog image generator
//  Copyright (C) 2013 Tom Lechner
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//



#include "catalogs.h"

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include </usr/include/GraphicsMagick/Magick++.h>
#include <sys/times.h>

#include <lax/strmanip.h>


#include <lax/laximages.h>


#include <lax/lists.cc>


using namespace std;
using namespace Laxkit;
using namespace LaxFiles;
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


#define DBG

//------------------------------- RenderContext -----------------------------------
/*! \class RenderContext
 * Info about how, what, and where to render.
 */



RenderContext::RenderContext()
{

	projectfile=newstr("save.stars");
	filename   =newstr("stars.tif");
	width      =8192;
	height     =4096;
	galactic   =0;
	transparent=0;


	maxstarsize  =20; //pixels of biggest star
	bigthreshhold=6; //magnitude < this get big treatment
	maximum_magnitude=0;
	minimum_magnitude=18;

	min_asc=0;
	max_asc=360;
	min_dec=-90;
	max_dec=90;

	usehalo  =5;
	halo     =NULL;
	halowidth=0;
	halotype=0;

	data=NULL;
	catalog=NULL;

	 //curve objects:
	ramp        =new CurveInfo("Halo Ramp", "radius",0,1,  "a",0,255);
	blowout     =new CurveInfo("Halo Blowout", "opacity",0,1,  "color",0,1);
	index_r     =new CurveInfo("Red", "B-V",-0.5,2,  "r",0,255);
	index_g     =new CurveInfo("Green", "B-V",-0.5,2,  "g",0,255);
	index_b     =new CurveInfo("Blue", "B-V",-0.5,2, "b",0,255);
	bigscale    =new CurveInfo("Big star scale", "magnitude",bigthreshhold,-1, "px",1,50);
	pointopacity=new CurveInfo("Point Opacity", "mag",minimum_magnitude,bigthreshhold, "a",0,1);

	Reset(~0);
}

RenderContext::~RenderContext()
{
	if (projectfile) delete[] projectfile;
	if (filename) delete[] filename;

	ramp        ->dec_count();
	blowout     ->dec_count();
	index_r     ->dec_count();
	index_g     ->dec_count();
	index_b     ->dec_count();
	bigscale    ->dec_count();
	pointopacity->dec_count();
}

void RenderContext::Reset(int which)
{

	if (which&POINTOPACITY) pointopacity->Reset();
	if (which&BIGSCALE) {
		bigscale->Reset();
		bigscale->MovePoint(1, 1,20);
	}

	if (which&RAMP) ramp->Reset();
	if (which&BLOWOUT) blowout->Reset();

	if (which&INDEX_R) {
		index_r->Reset();
		index_r->AddPoint(.25,230);
		index_r->AddPoint(.4,255);
	}

	if (which&INDEX_G) {
		index_g->Reset();
		index_g->MovePoint(0, -.5,170);
		index_g->MovePoint(1, 2,0);
		index_g->AddPoint(.4,255);
		index_g->AddPoint(1.7,200);
	}

	if (which&INDEX_B) {
		index_b->Reset();
		index_b->MovePoint(0, -.5,255);
		index_b->MovePoint(1, 2,0);
		index_b->AddPoint(.5,255);
		index_b->AddPoint(1.7,150);
	}
}

void RenderContext::InstallHaloImage(const char *file)
{
	cout <<"Installing halo file: "<<file<<endl;
	//cerr <<" *** FAIL for imlib based when no x yet!! must rewrite"<<endl;
	return;

	Image image(file);
	if (image.rows()!=image.columns()) {
		cerr <<" Halo image must be same width and height. Ignoring."<<endl;
		return;
	}

	halowidth=image.rows();

	halotype=1;
	if (halo) delete[] halo;
	halo=new unsigned char[halowidth*halowidth*4];

	image.write(0,0, halowidth,halowidth,"BGRA",CharPixel, halo);


//	-------imlib based:--------------------
//	LaxImage *img=load_image(file); // *** FAIL for imlib based when no x
//	if (!img) return;
//	if (img->w()!=img->h()) {
//		img->dec_count();
//		cerr <<" Halo image must be same width and height. Ignoring."<<endl;
//	}
//
//	if (halo) delete[] halo;
//	halotype=1;
//	unsigned char *hfile=img->getImageBuffer();
//
//	halowidth=img->w();
//	halo=new unsigned char[halowidth*halowidth*4];
//	memcpy(halo, hfile, halowidth*halowidth*4);
//
//	img->doneWithBuffer(hfile);
//	img->dec_count();


//  ---debug write out:
//	const char *saveto="halo.png";
//	if (saveto) {
//		Image haloimage;
//
//		int usecolor=1;
//		int w=halowidth;
//		if (usecolor) haloimage.read(w,w,"BGRA",CharPixel,halo);
//		else haloimage.read(w,w,"I",CharPixel,halo);
//
//		haloimage.magick("PNG");
//		haloimage.write(saveto);
//	}

}

void RenderContext::RefreshStats()
{
	stats.Zero();
	int curtime=times(NULL);

	for (int c=0; c<catalogs.n; c++) {
		catalog=catalogs.e[c];
		if (!catalog->visible) continue;

		catalogs.e[c]->RefreshStats(this,1);

		stats.numstars+=catalogs.e[c]->stats.numstars;
		stats.numgalaxies+=catalogs.e[c]->stats.numgalaxies;
		stats.numnebulae+=catalogs.e[c]->stats.numnebulae;
		stats.numother+=catalogs.e[c]->stats.numother;

		for (int c2=0; c2<50; c2++) stats.mags[c2]+=catalogs.e[c]->stats.mags[c2];

		if (catalogs.e[c]->stats.magmin<stats.magmin) stats.magmin=catalogs.e[c]->stats.magmin;
		if (catalogs.e[c]->stats.magmax>stats.magmax) stats.magmax=catalogs.e[c]->stats.magmax;
	}

	stats.rendertime=times(NULL)-curtime;
}



//! --- Actual Render and Export Function ---
int RenderContext::Render()
{
	DBG cout <<"Render with config:"<<endl;
	DBG dump_out(stdout,0,0,NULL);

	cout <<"Render to: "<<filename<<endl;
	cout <<"Width:     "<<width<<endl;
	cout <<"Height:    "<<height<<endl;

	//InstallHaloImage("smiley.png"); // ***

	 //create halo image if necessary
	if (!halo) {
		halowidth=maxstarsize*usehalo;
		halo=new unsigned char[halowidth*halowidth];
		halotype=0;
		CreateStockHalo(halowidth, usehalo, halo, "g",
						ramp,blowout,NULL);
						//ramp,blowout,"halo.png");
	}

	 //-----allocate star image data
	if (!data) {
		data=new unsigned char[(long)width*height*4];
	}

	int stride=width*4;
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



	int mags[50];
	int magmin=1000, magmax=-1000;
	memset(mags,0,50*sizeof(int));
	int numstars=0;
	int numgalaxies=0;

	 //--------- Render the stack of catalogs ---------------------
	for (int c=0; c<catalogs.n; c++) {
		catalog=catalogs.e[c];
		if (!catalog->visible) continue;
		catalog->Render(this);
	}



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


	stars->compressType(LZWCompression);
	stars->read(width,height,"BGRA",CharPixel,data);
    stars->magick("TIFF");


	cout <<"Writing to "<<filename<<"..."<<endl;
	stars->write(filename);


	cout << "Done!"<<endl;

	return 0;
}


void RenderContext::dump_out(FILE *f,int indent,int what,anObject *context)
{
    Attribute *att=dump_out_atts(NULL,0,context);
    att->dump_out(f,indent);
    delete att;
}

LaxFiles::Attribute *RenderContext::dump_out_atts(LaxFiles::Attribute *att,int what,anObject *context)
{
    if (!att) att=new Attribute("RenderContext",NULL);

    if (what==-1) {
        // *** format description
        return att;
    }

    if (filename) att->push("renderto",filename);
	att->push("width", width,-1);
	att->push("height",height,-1);

	//att->push("principal_galaxy_catalog","***");
	//att->push("tycho_2_star_catalog","***");

	att->push("galactic",galactic?"yes":"no");
	att->push("transparent",transparent?"yes":"no");
	
	att->push("bigthreshhold",bigthreshhold,-1);
	att->push("minimum_magnitude",minimum_magnitude,-1);
	att->push("maximum_magnitude",maximum_magnitude,-1);

	att->push("halosize",usehalo,-1);

	att->push("min_asc", min_asc,-1);
	att->push("max_asc", max_asc,-1);
	att->push("min_dec", min_dec,-1);
	att->push("max_dec", max_dec,-1);

	Attribute *aa;
	aa=new Attribute("ramp",NULL);
	ramp->dump_out_atts(aa,what,context);
	att->push(aa,-1);

	aa=new Attribute("blowout",NULL);
	blowout->dump_out_atts(aa,what,context);
	att->push(aa,-1);

	aa=new Attribute("index_r",NULL);
	index_r->dump_out_atts(aa,what,context);
	att->push(aa,-1);

	aa=new Attribute("index_g",NULL);
	index_g->dump_out_atts(aa,what,context);
	att->push(aa,-1);

	aa=new Attribute("index_b",NULL);
	index_b->dump_out_atts(aa,what,context);
	att->push(aa,-1);

	aa=new Attribute("bigscale",NULL);
	bigscale->dump_out_atts(aa,what,context);
	att->push(aa,-1);

	aa=new Attribute("pointopacity",NULL);
	pointopacity->dump_out_atts(aa,what,context);
	att->push(aa,-1);




    //char ss[50], *str=NULL;
	Catalog *cat;
    for (int c=0; c<catalogs.n; c++) {
		cat=catalogs.e[c];
		aa=new Attribute("catalog",cat->TypeName());
		cat->dump_out_atts(aa,what,context);
		att->push(aa,-1);
    }


    return att;
}

void RenderContext::dump_in_atts(LaxFiles::Attribute *att,int flag,anObject *context)
{
	const char *pgc_file_119=NULL;
	const char *pgc_file_237=NULL;
	const char *tycho_file=NULL;

    char *name,*value;
    for (int c=0; c<att->attributes.n; c++) {
        name= att->attributes.e[c]->name;
        value=att->attributes.e[c]->value;

        if (!strcmp(name,"filename")) {
            makestr(filename,value);

		} else if (!strcmp(name,"width")) {
			LongAttribute(value,&width,NULL);

		} else if (!strcmp(name,"height")) {
			LongAttribute(value,&height,NULL);

		} else if (!strcmp(name,"principal_galaxy_catalog_237")) {
			pgc_file_237=value;

		} else if (!strcmp(name,"principal_galaxy_catalog_119")) {
			pgc_file_119=value;

		} else if (!strcmp(name,"tycho_2_star_catalog")) {
			tycho_file=value;

		} else if (!strcmp(name,"galactic")) {
			galactic=BooleanAttribute(value);

		} else if (!strcmp(name,"transparent")) {
			transparent=BooleanAttribute(value);
	
		} else if (!strcmp(name,"bigthreshhold")) {
			DoubleAttribute(value,&bigthreshhold,NULL);

		} else if (!strcmp(name,"minimum_magnitude")) {
			DoubleAttribute(value,&minimum_magnitude,NULL);

		} else if (!strcmp(name,"maximum_magnitude")) {
			DoubleAttribute(value,&maximum_magnitude,NULL);

		} else if (!strcmp(name,"halosize")) {
			DoubleAttribute(value,&usehalo,NULL);

		} else if (!strcmp(name,"min_asc")) {
			DoubleAttribute(value,&min_asc,NULL);

		} else if (!strcmp(name,"max_asc")) {
			DoubleAttribute(value,&max_asc,NULL);

		} else if (!strcmp(name,"min_dec")) {
			DoubleAttribute(value,&min_dec,NULL);

		} else if (!strcmp(name,"max_dec")) {
			DoubleAttribute(value,&max_dec,NULL);



		} else if (!strcmp(name,"ramp")) {
			ramp->dump_in_atts(att->attributes.e[c], flag,context);

		} else if (!strcmp(name,"blowout")) {
			blowout->dump_in_atts(att->attributes.e[c], flag,context);

		} else if (!strcmp(name,"index_r")) {
			index_r->dump_in_atts(att->attributes.e[c], flag,context);

		} else if (!strcmp(name,"index_g")) {
			index_g->dump_in_atts(att->attributes.e[c], flag,context);

		} else if (!strcmp(name,"index_b")) {
			index_b->dump_in_atts(att->attributes.e[c], flag,context);

		} else if (!strcmp(name,"bigscale")) {
			bigscale->dump_in_atts(att->attributes.e[c], flag,context);

		} else if (!strcmp(name,"pointopacity")) {
			pointopacity->dump_in_atts(att->attributes.e[c], flag,context);

		} else if (!strcmp(name,"catalog")) {
			Catalog *cat=NULL;
			if (!strcmp(value,"RandomMemory")) cat=new RandomCatalog(NULL,1,0);
			else cat=new Catalog(NULL,NULL,Unknown);

			cat->dump_in_atts(att->attributes.e[c], flag,context);
			catalogs.push(cat,1);
		}
	}

	bigscale    ->RefreshLookup(256, 1,bigscale->ymax);
	pointopacity->RefreshLookup(256, 0,255);
	ramp        ->RefreshLookup(256, 0,255);
	blowout     ->RefreshLookup(256, 0,255);
	index_r     ->RefreshLookup(256, 0,255);
	index_g     ->RefreshLookup(256, 0,255);
	index_b     ->RefreshLookup(256, 0,255);



	if (pgc_file_119) catalogs.push(new Catalog("Principal Galaxies Catalog (119)", pgc_file_119, PrincipalGalaxy119),1);
	if (pgc_file_237) catalogs.push(new Catalog("Principal Galaxy Catalog (237)", pgc_file_237, PrincipalGalaxy237),1);
	if (tycho_file)   catalogs.push(new Catalog("Tycho 2 Star Catalog",     tycho_file,   Tycho2),            1);
}

 

//------------------------------- StarPoint -----------------------------------

void StarPoint::Set(double nindex, double nmag, double nasc, double ndec, int ntype)
{
	type=ntype;
	index=nindex;
	mag=nmag;
	asc=nasc;
	dec=ndec;
}


//------------------------------- Random Catalog for previewing mainly -----------------------------------


/*! \class RandomCatalog
 * Collection of stars in memory for previewing.
 */
RandomCatalog::RandomCatalog(const char *nname, int num, int spherical)
  : Catalog(nname,NULL,RandomMemory)
{
	isspherical=spherical;

	if (num<=0) num=1000;

	if (spherical==2 || spherical==3) RepopulateFakeMilkyWay(num,spherical==2?0:1);
	else Repopulate(num,isspherical);
}

//! Convert cartesian space to spherical surface. (radians)
void rect_to_sphere(double x,double y,double z, double *asc,double *dec)
{
	*asc=atan2(y,x);
	*dec=atan(z/sqrt(x*x+y*y));
}

int RandomCatalog::RepopulateFakeMilkyWay(int num,int galactic)
{
	isspherical=1;
	num_cat_points=0;

	 //make fake milky way
	double asc,dec;
	flatpoint p;
	double dd=1.5;

	for (int c=0; c<num; c++) {
		if (c>=points.n) {
			points.push(new StarPoint(0,0,0,0,0));
			num_cat_points=c+1;
		}

		//normalized:
		asc=drand48(); //this creates bunching at the poles
		dec=tan(drand48()*2*dd-dd)/tan(dd); //-1..1
		dec=dec*.5+.5;

		//if (!galactic) {
			//p=Gal2Eq(asc*2*M_PI, dec*M_PI-M_PI/2);
			//asc=p.x/360;
			//dec=p.y/180+.5;
		//}

		//spherical:
		//asc=drand48()*360;
		//dec=tan(drand48()*2*dd-dd)/tan(dd)*180 - 90;

		points[c]->asc=asc;
		points[c]->dec=dec;
		points[c]->index=drand48()*(2.5)-.5;
		points[c]->mag  =drand48()*(15);
	}

	return num_cat_points;
}

int RandomCatalog::Repopulate(int num, int spherical)
{
	isspherical=spherical;


	double x,y,z;
	num_cat_points=0;
	for (int c=0; c<num; c++) {
		if (c>=points.n) {
			points.push(new StarPoint(0,0,0,0,0));
			num_cat_points=c+1;
		}

		points[c]->index=drand48()*(2.5)-.5;
		points[c]->mag  =drand48()*(15);
		
		 //create random point in unit cube, map to sphere..
		 //creates even distribution compared to just random asc/dec
		if (spherical) {
			x=drand48()*2-1;
			y=drand48()*2-1;
			z=drand48()*2-1;

			if (x*x+y*y+z*z>1) { c--; continue; } //make a smooth spherical distribution

			if (x==0 && y==0 && z==0) x=1;
			rect_to_sphere(x,y,z, &points[c]->asc, &points[c]->dec);
			//makes asc [-180..180], dec [-90..90] but in radians

			points[c]->asc/=2*M_PI; //convert to unit range
			points[c]->asc+=.5;
			points[c]->dec/=M_PI; //convert to unit range
			points[c]->dec+=.5;


		} else {
			 //for simple rectilinear:
			points[c]->asc=drand48();
			points[c]->dec=drand48();
		}
	}

	return num_cat_points;
}

int RandomCatalog::Render(RenderContext *context)
{
    double Ra;   //right ascension (longitude)
    double Dec;  //declination (latitude)
    double vmag; //visual mag
    double index;
	double bmag;


	for (int c=0; c<num_cat_points; c++) {

		vmag  = points.e[c]->mag;
		index = points.e[c]->index;
		bmag  = index+vmag;
		Ra    = points.e[c]->asc*360;
		Dec   = points.e[c]->dec*180-90;

		if (vmag>context->minimum_magnitude || vmag<context->maximum_magnitude) continue;

		drawStar(context, Ra, Dec, vmag, bmag);
	}

	return 0;
}

int RandomCatalog::Render(RenderContext *context, unsigned char *data,int ww,int hh)
{
	int curtime=times(NULL);

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

	for (int c=0; c<num_cat_points; c++) {

		vmag  = points.e[c]->mag;
		index = points.e[c]->index;
		Ra    = points.e[c]->asc;
		Dec   = points.e[c]->dec;

		if (vmag>context->minimum_magnitude || vmag<context->maximum_magnitude) continue;

		drawStarSimple(context, Ra, Dec, index, vmag);
	}


	 //restore actual surface
	context->data  =olddata;
	context->width =oldw;
	context->height=oldh;  

	stats.rendertime=times(NULL)-curtime;

	return 0;
}

int RandomCatalog::RefreshStats(RenderContext *context, int buildpoints)
{
	int curtime=times(NULL);

    //double Ra;   //right ascension (longitude)
    //double Dec;  //declination (latitude)
    double vmag; //visual mag
    //double index;//blue mag

	int numstars=0;

	double magmin=-5000;
	double magmax=5000;
	int mags[50];
	memset(mags,0,50*sizeof(int));

	for (int c=0; c<num_cat_points; c++) {

		vmag  = points.e[c]->mag;
		//index = points.e[c]->index;
		//Ra    = points.e[c]->asc;
		//Dec   = points.e[c]->dec;

		if (vmag>context->minimum_magnitude || vmag<context->maximum_magnitude) continue;

		if (vmag<magmin) magmin=vmag;
		if (vmag>magmax) magmax=vmag;
		mags[int(vmag>=0?vmag:0)]++;

		numstars++;

	}

	stats.Zero();
	stats.numstars=numstars;
	memcpy(stats.mags, mags, 50*sizeof(int));
	stats.magmin=magmin;
	stats.magmax=magmax;

	stats.rendertime=times(NULL)-curtime;

	return 0;
}


//------------------------------- CatalogStats -----------------------------------

/*! \class CatalogStats
 */

CatalogStats::CatalogStats()
{
	numstars=-1;
	numgalaxies=-1;
	numnebulae=-1;
	numother=-1;

	memset(mags,0,50*sizeof(int));
	zeromag=0; //what is the actual magnitude of mags[0], if mag[0] is negative;

	magmin=1000;
	magmax=-1000;

	rendertime=0;
}

//! Zero out all stats. Initially, all is -1 for unknown.
void CatalogStats::Zero()
{
	numstars   =0;
	numgalaxies=0;
	numnebulae =0;
	numother   =0;

	memset(mags,0,50*sizeof(int));
	zeromag=0; //what is the actual magnitude of mags[0], if mag[0] is negative;

	magmin=1000;
	magmax=-1000;
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
	//magnitude_distribution=NULL;
	//nummags=0;
	minimum_magnitude=maximum_magnitude=0;

	min_mag_cutoff=20;
	max_mag_cutoff=0;

	 //StarPoint takes 36+pointer=44 bytes per star
	 //so 4m == 44 megs. This is an arbitrary limit, just depends on user memory
	max_points_in_memory=4000000;
	num_cat_points=0;

	visible=1;
	autoadded=0;

	points.Delta(20000);
}

Catalog::~Catalog()
{
	if (name) delete[] name;
	if (filename) delete[] filename;
}

const char *Catalog::TypeName()
{
	if (type==Tycho2)          return "Tycho2";
	if (type==PrincipalGalaxy119) return "PrincipalGalaxy119";
	if (type==PrincipalGalaxy237) return "PrincipalGalaxy237";
	if (type==RandomMemory)    return "RandomMemory";
	if (type==CustomStars)     return "CustomStars";
	if (type==CustomGalaxy)    return "CustomGalaxy";
	return "Unknown";
}

void Catalog::dump_out(FILE *f,int indent,int what,anObject *context)
{
    Attribute *att=dump_out_atts(NULL,0,context);
    att->dump_out(f,indent);
    delete att;
}

LaxFiles::Attribute *Catalog::dump_out_atts(LaxFiles::Attribute *att,int what,anObject *context)
{
    if (!att) att=new Attribute("RenderContext",NULL);

    if (what==-1) {
        // *** format description
        return att;
    }

    if (name) att->push("name",name);
	att->push("type",TypeName());
    if (filename) att->push("file",filename);

	att->push("minimum_magnitude",minimum_magnitude,-1);
	att->push("maximum_magnitude",maximum_magnitude,-1);
	att->push("max_mag_cutoff",max_mag_cutoff,-1);
	att->push("min_mag_cutoff",min_mag_cutoff,-1);

	//att->push("min_asc", min_asc,-1);
	//att->push("max_asc", max_asc,-1);
	//att->push("min_dec", min_dec,-1);
	//att->push("max_dec", max_dec,-1);

	//Attribute *aa;
	//aa=new Attribute("ramp",NULL);
	//ramp->dump_out_atts(aa,what,context);
	//att->push(aa,-1);

    return att;
}

void Catalog::dump_in_atts(LaxFiles::Attribute *att,int flag,anObject *context)
{
    char *aname,*value;
    for (int c=0; c<att->attributes.n; c++) {
        aname=att->attributes.e[c]->name;
        value=att->attributes.e[c]->value;

        if (!strcmp(aname,"file")) {
            makestr(filename,value);

		} else if (!strcmp(aname,"name")) {
            makestr(name,value);

		} else if (!strcmp(aname,"type")) {
			if (!strcmp(value,"Tycho2")) type=Tycho2;
			else if (!strcmp(value,"PrincipalGalaxy119")) type=PrincipalGalaxy119;
			else if (!strcmp(value,"PrincipalGalaxy237")) type=PrincipalGalaxy237;
			//other types are subclasses...?

		} else if (!strcmp(aname,"minimum_magnitude")) {
			DoubleAttribute(value,&minimum_magnitude,NULL);

		} else if (!strcmp(aname,"maximum_magnitude")) {
			DoubleAttribute(value,&maximum_magnitude,NULL);

		} else if (!strcmp(aname,"min_mag_cutoff")) {
			DoubleAttribute(value,&min_mag_cutoff,NULL);

		} else if (!strcmp(aname,"max_mag_cutoff")) {
			DoubleAttribute(value,&max_mag_cutoff,NULL);

		//} else if (!strcmp(name,"pointopacity")) {
		//	pointopacity->dump_in_atts(att->attributes.e[c], flag,context);
		
		}

	}

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
	if (type==CustomGalaxy)    return Process_Galaxy(context,0,0);
	if (type==PrincipalGalaxy119) return Process_PGC(context,0,0, 119);
	if (type==PrincipalGalaxy237) return Process_PGC(context,0,0, 237);
	if (type==Tycho2)          return Process_Tycho(context,0,0);

	return 1;
}

int Catalog::RefreshStats(RenderContext *context, int buildpoints)
{
	if (type==CustomGalaxy)    return Process_Galaxy(context,1, buildpoints);
	if (type==PrincipalGalaxy119) return Process_PGC(context,1, buildpoints, 119);
	if (type==PrincipalGalaxy237) return Process_PGC(context,1, buildpoints, 237);
	if (type==Tycho2)          return Process_Tycho(context,1, buildpoints);

	return 0;
}




//------------------------- Custom Galaxy processing ---------------------------
int Catalog::Process_Galaxy(RenderContext *rr, int statsonly, int buildpoints)
{
	const char *file=filename;

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

	int curtime=times(NULL);

	cout <<"Scanning Galaxy file: "<<file<<endl;

	int numstars=0;
	if (file) {
		char *line=NULL;
		size_t n=0;
		ssize_t c;
		FILE *f=fopen(file,"r");
		if (!f) {
			cerr <<" --Fail!-- Could not open Custom Galaxy file: "<<file<<endl;
			return 1;
		}
		do {

			c=getline(&line,&n,f);
			if (c<0) break;
			if (c==0) continue;

			 //the catalog lines are delimited by '|' characters, but also arranged on strict byte widths
			Ra   = strtod(line,NULL);          // RA
			Dec  = strtod(line+13,NULL);         //DEC
			vmag = strtod(line+26,NULL); //overall visual magnitude, bytes 60-63
			vmag-= 4;
			bmag = vmag+.2; //makes index refer to a slightly blue object

			if (Ra<0) Ra+=360;
		 
			if (Ra<rr->min_asc || Ra>rr->max_asc || Dec<rr->min_dec || Dec>rr->max_dec) continue;

			if (buildpoints && numstars<max_points_in_memory) {
				if (numstars<points.n) points.e[numstars]->Set(bmag-vmag, vmag, Ra, Dec, 0);
				else points.push(new StarPoint(bmag-vmag, vmag, Ra, Dec, 0),1);
				num_cat_points=numstars;
			}

			//cerr <<" a,d,m: "<<Ra<<"  "<<Dec<<"  "<<vmag<<endl;

			if (vmag>rr->minimum_magnitude || vmag<rr->maximum_magnitude) continue;
			if (Ra==0 && Dec==0) continue;

			if (vmag<magmin) magmin=vmag;
			if (vmag>magmax) magmax=vmag;
			mags[int(vmag>=0?vmag:0)]++;


			if (!statsonly) drawStar(rr, Ra, Dec, vmag, bmag);

			numstars++;
			if (numstars%10000==0) cout <<"+\n";

		} while (!feof(f));
		if (line) free(line);

		fclose(f);
	}

	stats.Zero();
	stats.numstars=numstars;
	memcpy(stats.mags, mags, 50*sizeof(int));
	stats.magmin=magmin;
	stats.magmax=magmax;

	stats.rendertime=times(NULL)-curtime;

	return 0;
}



//------------------------- Principal Galaxy Catalog processing ---------------------------
int Catalog::Process_PGC(RenderContext *rr, int statsonly, int buildpoints, int revision)
{
	const char *pgc_file=filename;
	if (!pgc_file) return 1;


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

	int curtime=times(NULL);

	//default version 119 values:
	int byte_Ra=6;
	int byte_Dec=14;
	int byte_mag=59;
	//int byte_type=39;
	//int byte_major=43;
	//int byte_minor=51;
	//int byte_angle=73;

	if (revision==237) {
		byte_Ra=12;
		byte_Dec=20;
		//byte_type=31;
		//byte_major=36; //label: logD25, units: 0.1arcmin
		//byte_minor=50; //label: logR25, ==log(majoraxis/minoraxis)
		//byte_angle=63;//

		byte_mag=-1; //*****major axis has: Apparent diameter in log scale (D in 0.1arcmin) converted to the
				     //                     RC3 system at the isophote 25B-mag/arcsec^2^ (section 3);
					 //an isophote is a contour line of the same brightness
	} 





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
			Ra   = hms(line+byte_Ra);          // RA hours-min-sec,        bytes 7-14
			Dec  = dms(line+byte_Dec);         //DEC degrees-min-sec,      bytes 15-21
			vmag = strtod(line+byte_mag,NULL); //overall visual magnitude, bytes 60-63
			bmag = vmag+.2; //makes index refer to a slightly blue object

			// *** type of galaxy, bytes 40-43
			// *** major axis, arcmin, bytes 44-49
			// *** minor axis, arcmin, bytes 52-56
			// *** Position Angle from North eastward, bytes 74-76
		 
			if (vmag>rr->minimum_magnitude || vmag<rr->maximum_magnitude) continue;
			if (Ra==0 && Dec==0) continue;

			if (Ra<rr->min_asc || Ra>rr->max_asc || Dec<rr->min_dec || Dec>rr->max_dec) continue;

			if (vmag<magmin) magmin=vmag;
			if (vmag>magmax) magmax=vmag;
			mags[int(vmag>=0?vmag:0)]++;

			if (buildpoints && numgalaxies<max_points_in_memory) {
				if (numgalaxies<points.n) points.e[numgalaxies]->Set(bmag-vmag, vmag, Ra, Dec, 0);
				else points.push(new StarPoint(bmag-vmag, vmag, Ra, Dec, 0),1);
				num_cat_points=numgalaxies;
			}


			if (!statsonly) drawStar(rr, Ra, Dec, vmag, bmag);

			numgalaxies++;
			if (numgalaxies%10000==0) cout <<"+\n";

		} while (!feof(f));
		if (line) free(line);

		fclose(f);
	}
	cout <<"Rendered "<<numgalaxies<<" from "<<pgc_file<<endl;

	stats.Zero();
	stats.numgalaxies=numgalaxies;
	memcpy(stats.mags, mags, 50*sizeof(int));
	stats.magmin=magmin;
	stats.magmax=magmax;

	stats.rendertime=times(NULL)-curtime;

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

int Catalog::Process_Tycho(RenderContext *rr, int statsonly, int buildpoints)
{
	const char *tycho_file=filename;



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

	int curtime=times(NULL);


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

			if (Ra<rr->min_asc || Ra>rr->max_asc || Dec<rr->min_dec || Dec>rr->max_dec) continue;

			if (vmag<magmin) magmin=vmag;
			if (vmag>magmax) magmax=vmag;
			mags[int(vmag-magmin)]++;

			 //find stats about colors:
			//colorindex=bmag-vmag;
			//if (colorindex<colorindex_min) colorindex_min=colorindex;
			//if (colorindex>colorindex_max) colorindex_max=colorindex;

			if (buildpoints && numstars<max_points_in_memory) {
				if (numstars<points.n) points.e[numstars]->Set(bmag-vmag, vmag, Ra, Dec, 0);
				else points.push(new StarPoint(bmag-vmag, vmag, Ra, Dec, 0),1);
				num_cat_points=numstars;
			}


			if (!statsonly) drawStar(rr, Ra, Dec, vmag, bmag);

			numstars++;
			if (numstars%100000==0) cout <<".\n";

		} while (!feof(f));
		if (line) free(line);

		fclose(f);
	}

	stats.Zero();
	stats.numstars=numstars;
	memcpy(stats.mags, mags, 50*sizeof(int));
	stats.magmin=magmin;
	stats.magmax=magmax;

	stats.rendertime=times(NULL)-curtime;

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
	int b=context->data[i+0];
	int g=context->data[i+1];
	int r=context->data[i+2];
	int a=context->data[i+3];

	double aaa=color.alphaf();
	int aa=aaa*255;
	int rr=aaa*color.red();
	int gg=aaa*color.green();
	int bb=aaa*color.blue();

	r+=rr; if (r>255) r=255;
	g+=gg; if (g>255) g=255;
	b+=bb; if (b>255) b=255;
	a+=aa; if (a>255) a=255;

	context->data[i  ]=b;
	context->data[i+1]=g;
	context->data[i+2]=r;
	context->data[i+3]=a;

	//cerr <<"blend rgba  "<<r<<' '<<g<<' '<<b<<' '<<a<<endl;
}

/*! New color= (oldcolor)*(1-a) + newcolor*a
 */
void alphaOverPixel(RenderContext *context, int x,int y, StarColor &color)
{
	if (x<0 || x>=context->width || y<0 || y>=context->height) return;

	int i=y*context->width*4+x*4;
	int b=context->data[i+0];
	int g=context->data[i+1];
	int r=context->data[i+2];
	int a=context->data[i+3];

	double aaa=color.alphaf();
	int aa=aaa*255;
	int rr=aaa*color.red();
	int gg=aaa*color.green();
	int bb=aaa*color.blue();

	r=r*(1-aaa) + rr*aaa;
	g=g*(1-aaa) + gg*aaa;
	b=b*(1-aaa) + bb*aaa;
	a+=aa; if (a>255) a=255;

	context->data[i  ]=b;
	context->data[i+1]=g;
	context->data[i+2]=r;
	context->data[i+3]=a;

	//cerr <<"blend rgba  "<<r<<' '<<g<<' '<<b<<' '<<a<<endl;
}

void dataEllipse(RenderContext *context, int xp,int yp, double xr,double yr, StarColor &color)
{
	double xspan;

	if (context->usehalo) {
		 //may halo reference to data
		 //rendered on in a very naive way, just rectangular copy, no additional span correction

		StarColor col;
		int i;
		int aa; //blowout mapped alpha
		double a,r,g,b;
		int sx,sy; //offset into destintion rectangle
		int hx,hy; //coord in halo sample image

		for (int y=yp-yr; y<yp+yr; y++) {
		  sy=y-(yp-yr);
		  for (int x=xp-xr; x<xp+xr; x++) {
		    sx=x-(xp-xr);

			hy=sy/(2.*yr)*context->halowidth;
			hx=sx/(2.*xr)*context->halowidth;
			i=hy*context->halowidth + hx;

			if (context->halotype!=0) {
				 //if custom color halo image, map image in alpha over
				i=i*4;
				b =context->halo[i++];
				g =context->halo[i++];
				r =context->halo[i++];
				a =context->halo[i++];
				col.alpha(a);

				col.red  (r);
				col.green(g);
				col.blue (b);

				alphaOverPixel(context, x,y, col);
			} else {
				 //normal additive star blending
				r=color.redf();
				g=color.greenf();
				b=color.bluef();

				a=context->halo[i]/255.;//so this is the halo alpha 0..1

				col.alpha(a*color.alpha());
				aa=context->blowout->lookup[int(a*255)];

				r=r*(255-aa) + (aa);
				g=g*(255-aa) + (aa);
				b=b*(255-aa) + (aa);
			
				col.red  (r);
				col.green(g);
				col.blue (b);

				blendPixel(context, x,y, col);
			}

		  }
		  //cecontext <<endl;
		}
		

	} else {
		 //draw (alas non-antialiased) pixels in a circle
		for (int y=yp-yr; y<yp+yr; y++) {
		  xspan=xr*sqrt(1-(y-yp)*(y-yp)/yr/yr);
		  for (int x=xp-xspan; x<xp+xspan; x++) {
			blendPixel(context, x,y, color);
		  }
		}
	}
}


void point(RenderContext *context, double x,double y, StarColor &color, double span)
{
	if (span<=1) {
		 //near enough to equator that it is just a single pixel
		blendPixel(context, x,y,color);

	} else if (x-span/2<0) {
		 //draw partial, wraps around 
		for (int c=0; c<x+span/2; c++) blendPixel(context, c,y,color); //draw line 0 to x+span/2
		for (int c=x+context->width-span/2; c<context->width; c++) blendPixel(context, c,y,color); //draw line x+width-span/2, to width

	} else if (x+span/2>context->width) {
		 //wraps around to beginning
		for (int c=0; c<x-context->width+span/2; c++) blendPixel(context, c,y,color); //line 0 to x-width+span/2
		for (int c=x-span/2; c<context->width; c++) blendPixel(context, c,y,color); //line x-span/2 to width

	} else {
		 //fits without wraparound
		for (int c=x-span/2; c<x+span/2; c++) blendPixel(context, c,y,color); //line x-span/2 to width
	}
}


/*! span is how much horizontally to stretch out one pixel.
 * Vertical is not stretched. Note this is not accurate.. ellipses do not actually project into ellipses.
 */
void ellipse(RenderContext *rr, double x,double y, double xs,double ys, StarColor &color, double span)
{

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
 * ra and dec are assumed to be in range 0..1, and mapped to context->width and height.
 */
void drawStarSimple(RenderContext *context, double ra, double dec, double index, double vmag)
{
  double x, y;
  StarColor color;
  indexToRgb(context, index, vmag, color);  // Get a color for the star
  
  // map to screen coordinates
  x = map(ra,  0, 1,     0,   context->width);
  y = map(dec, 0, 1, context->height, 0);

  //cout <<"x:"<<x<<"  y:"<<y<<endl;
  
  
   //figure out how much the star is stretched out in final image
   //simple rectilinear has 1:1 always
  double span=1;


  // For bright stars draw a cicle, otherwise just a pixel
  if (vmag < context->bigthreshhold)
  {
    //double s = map(vmag, context->bigthreshhold, -1, 2, context->maxstarsize);

    double s = map(vmag, context->bigthreshhold, -1, 0, 255);
	int ss=context->bigscale->lookup[(int)constrain(s, 0,255)];
	//cerr <<" map size mag:"<<vmag<<"  i:"<<(int)constrain(s, 0,255)<<"  size:"<<ss<<endl;

	color.alpha(color.alpha()/4+192); // *** this needs its own ramp
	ellipse(context, x, y, ss, ss,  color, span);
  }
  else
  {
	int aa=context->pointopacity->lookup[color.alpha()];
	color.alpha(aa);
	point(context, x, y, color, span);
  }
}



/**
 * This draws a pixel or small circle to the screen.
 * ra and dec are in range [0..360] and [-90..90], mapped to context->width and height.
 */
void drawStar(RenderContext *context, double ra, double dec, double vmag, double bmag)
{
  double x, y;
  StarColor color;
  double index = bmag-vmag;
  indexToRgb(context, index, vmag, color);  // Get a color for the star
  
  // map to screen coordinates
  x = map(ra, 0, 360, 0, context->width);
  y = map(dec, -90, 90, context->height, 0);

  //cout <<"x:"<<x<<"  y:"<<y<<endl;
  
  
  // if Galactic Coordinates..
  if (context->galactic) {
  	  //cout <<"ra:"<<radians(ra)<<"  dec:"<<radians(dec)<<endl;

	  flatvector galacticCoord = Eq2Gal(radians(ra), radians(dec));
	  x = galacticCoord.x;
	  
	  // Put the center of the Galaxy at the center of the image
	  if(x > 180)
		x -= 360;
	  
	  // map to screen coordinates
	  x = map(x, -180, 180, 0, context->width);
	  y = map(galacticCoord.y, -90, 90, context->height, 0);
  }
  
   //figure out how much the star is stretched out in final image
  double span=1;
  if (context->galactic) dec=180.*(y-context->height/2)/context->height;
  double declination_radians=dec/180.*M_PI;
  if (declination_radians==M_PI) span=context->width;
  else span=1/cos(declination_radians);
  if (span>context->width) span=context->width;


  // For bright stars draw a cicle, otherwise just a pixel
  if (vmag < context->bigthreshhold)
  {
    //double s = map(vmag, context->bigthreshhold, -1, 2, context->maxstarsize);

    double s = map(vmag, context->bigthreshhold, -1, 0, 255);
	int ss=context->bigscale->lookup[(int)constrain(s, 0,255)];
	//cerr <<" map size mag:"<<vmag<<"  i:"<<(int)constrain(s, 0,255)<<"  size:"<<ss<<endl;

	color.alpha(color.alpha()/4+192); // *** this needs its own ramp
	ellipse(context, x, y, ss, ss,  color, span);
  }
  else
  {
	int aa=context->pointopacity->lookup[color.alpha()];
	color.alpha(aa);
	point(context, x, y, color, span);
  }
}


/**
 * Regarding numbers for default index -> color:
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
void indexToRgb(RenderContext *context, double index, double vmag, StarColor &color)
{
  int i=map(index, -.5,2, 0,255);
  i=constrain(i, 0,255);

  int r = context->index_r->lookup[i];
  int g = context->index_g->lookup[i];
  int b = context->index_b->lookup[i];


   //map magnitude range of point starsto 0..1. Big stars map to alpha=255
  int bright = map(vmag,  context->minimum_magnitude, context->bigthreshhold+.0001,  0, 255);
  bright=constrain(bright, 0,255);
  

  color.red  (r);
  color.green(g);
  color.blue (b);
  color.alpha(bright);
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
double an          = radians(32.93192);   // Galactic long of asc node on equator
double ngpRa       = radians(192.85948);  // RA of North Galactic Pole
double ngpDec      = radians(27.12825);   // Dec of North Galactic Pole
double cos_ngpDec  = cos(ngpDec);
double sin_ngpDec  = sin(ngpDec);
static double SMALL = 1e-20;

/** radians in, degrees out.
 * Convert Equtorial Coordinates to Galactic Coordinates
 * Based on code from libastro. You are not expected to understand this.
 */
flatvector Eq2Gal(double ra, double dec)
{
  double sin_dec, cos_dec, a, cos_a, sin_a, b, square, c, d; 
  double lat_gal;
  double lon_gal;
  
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


double pole_ra = radians(192.859508);
double pole_dec = radians(27.128336);
double posangle = radians(122.932-90.0);

//! Radians inputs, decimal outputs.
flatvector Gal2Eq(double l, double b)
{
	 //North galactic pole (J2000) -- according to Wikipedia

	double ra = atan2( (cos(b)*cos(l-posangle)), (sin(b)*cos(pole_dec) - cos(b)*sin(pole_dec)*sin(l-posangle)) ) + pole_ra;
	double dec = asin( cos(b)*cos(pole_dec)*sin(l-posangle) + sin(b)*sin(pole_dec) );
	return flatpoint(degrees(ra), degrees(dec));
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
					 const char *format, Laxkit::CurveInfo *ramp, Laxkit::CurveInfo *blowout,
					 const char *saveto)
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


	if (saveto) {
		Image haloimage;

		if (usecolor) haloimage.read(w,w,"BGRA",CharPixel,halodata);
		else haloimage.read(w,w,"I",CharPixel,halodata);

		haloimage.magick("PNG");
		haloimage.write(saveto);
	}
}






