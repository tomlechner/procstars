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



#ifndef CATALOGS_H
#define CATALOGS_H



#include <lax/lists.h>
#include <lax/vectors.h>
#include <lax/curvewindow.h>



//------------------------------- Color -----------------------------------
class StarColor
{
  public:
	int r,g,b,a;
	int    red () { return r; }
	double redf() { return r/255.; }
	void   redf(double rr) { r=rr*255; }
	void   red (int rr)    { r=rr; }

	int    green () { return g; }
	double greenf() { return g/255.; }
	void   greenf(double gg) { g=gg*255; }
	void   green (int gg)    { g=gg; }

	int    blue () { return b; }
	double bluef() { return b/255.; }
	void   bluef(double bb) { b=bb*255; }
	void   blue (int bb)    { b=bb; }

	int    alpha () { return a; }
	double alphaf() { return a/255.; }
	void   alphaf(double aa) { a=aa*255; }
	void   alpha (int aa)    { a=aa; }
};

//------------------------------- RenderContext -----------------------------------
class Catalog;

enum RCurveType {
	RAMP        =(1<<0),
	BLOWOUT     =(1<<1),
	INDEX_R     =(1<<2),
	INDEX_G     =(1<<3),
	INDEX_B     =(1<<4),
	BIGSCALE    =(1<<5),
    POINTOPACITY=(1<<6)
};

//------------------------------- CatalogStats -----------------------------------

class CatalogStats
{
  public:
	int numstars;
	int numgalaxies;
	int numnebulae;
	int numother;

	int mags[50];
	int zeromag;
	double magmin;//by real number, not brightness
	double magmax;//by real number, not brightness

	double rendertime;
	
	CatalogStats();
	void Zero();
};

class RenderContext : public LaxFiles::DumpUtility
{
  public:
	char *projectfile; //save settings here
	char *filename; //final stars file
	long width;
	long height;
	unsigned char *data;

	int galactic;
	int transparent;
	double magnitude; //initial hint for magnitude cutoff
	double maxmagnitude;
	double maxstarsize; //pixels of biggest star
	double bigthreshhold; //magnitude < this get big treatment

	double usehalo;
	int halotype; //0 for bw, 1 for argb
	unsigned char *halo;
	int halowidth; //width in pixels of halo image

	Laxkit::CurveInfo *ramp;
	Laxkit::CurveInfo *blowout;
	Laxkit::CurveInfo *index_r;
	Laxkit::CurveInfo *index_g;
	Laxkit::CurveInfo *index_b;
	Laxkit::CurveInfo *bigscale;
	Laxkit::CurveInfo *pointopacity;

	double min_asc;
	double max_asc;
	double min_dec;
	double max_dec;


	double minimum_magnitude; //dimmest (numerically > max)
	double maximum_magnitude; //brightest

	Laxkit::PtrStack<Catalog> catalogs;
	Catalog *catalog; //current, set sequentially to each catalogs during Render()

	CatalogStats stats;
	virtual void RefreshStats();
	virtual int Render();

	RenderContext();
	virtual ~RenderContext();
	virtual void Reset(int which);
	virtual void InstallHaloImage(const char *file);

	virtual void dump_out(FILE *f,int indent,int what,Laxkit::anObject *context);
    virtual LaxFiles::Attribute *dump_out_atts(LaxFiles::Attribute *att,int what,Laxkit::anObject *context);
    virtual void dump_in_atts(LaxFiles::Attribute *att,int flag,Laxkit::anObject *context);
};



//------------------------------- Catalog -----------------------------------
enum CatalogTypes
{
	Unknown,
	Tycho2,
	PrincipalGalaxy119,
	PrincipalGalaxy237,
	RandomMemory,
	CustomStars,
	CustomGalaxy
};

class StarPoint
{
  public:
	int type;
	double index;
	double mag;
	double asc;
	double dec;

	StarPoint(double nindex, double nmag, double nasc, double ndec, int ntype=0)
		: type(ntype), index(nindex), mag(nmag), asc(nasc), dec(ndec)
	  {}
	void Set(double nindex, double nmag, double nasc, double ndec, int ntype);

 	 //galaxy info:
	//double majoraxis;
	//double minoraxis;
	//double tilt;
};

class Catalog
{
  protected:
  public:
	char *name;
	char *filename;
	CatalogTypes type;
	CatalogStats stats;
	
	int visible;
	int autoadded;

	 //stats
	double minimum_magnitude;
	double maximum_magnitude;

	 //subset
	double min_mag_cutoff;
	double max_mag_cutoff;

	int max_points_in_memory;
	int num_cat_points;
	Laxkit::PtrStack<StarPoint> points;


	Catalog(const char *nname, const char *nfile, CatalogTypes ntype);
	virtual ~Catalog();

	virtual const char *TypeName();
	virtual int Render(RenderContext *context);
	virtual int RefreshStats(RenderContext *context, int buildpoints);

	virtual int OpenCatalog();
	//virtual int GetLine(int &object_type, double &asc, double &dec, double &bmag, double &vmag) = 0;
	virtual int CloseCatalog();

	virtual void dump_out(FILE *f,int indent,int what,Laxkit::anObject *context);
    virtual LaxFiles::Attribute *dump_out_atts(LaxFiles::Attribute *att,int what,Laxkit::anObject *context);
    virtual void dump_in_atts(LaxFiles::Attribute *att,int flag,Laxkit::anObject *context);


	int Process_Galaxy(RenderContext *rr, int statsonly, int buildpoints);
	int Process_PGC(RenderContext *context, int statsonly, int buildpoints, int revision);
	int Process_Tycho(RenderContext *context, int statsonly, int buildpoints);
};

class RandomCatalog : public Catalog
{
  public:
	int isspherical;

	RandomCatalog(const char *nname, int num, int spherical);
	virtual int Render(RenderContext *context);
	virtual int Render(RenderContext *context, unsigned char *data,int ww,int hh);
	virtual int RefreshStats(RenderContext *context, int buildpoints);

	virtual int Repopulate(int num, int spherical);
	virtual int RepopulateFakeMilkyWay(int num,int galactic);
};

//------------------------------- Rendering Misc -----------------------------------

flatvector Eq2Gal(double ra, double dec);
flatvector Gal2Eq(double l, double b);
void indexToRgb(RenderContext *rr, double index, double vmag, StarColor &color);
void drawStar(RenderContext *rr, double ra, double dec, double vmag, double bmag);
void drawStarSimple(RenderContext *context, double ra, double dec, double index, double vmag);
void CreateStockHalo(int w,double halosize, unsigned char *halo, const char *format,
					 Laxkit::CurveInfo *ramp, Laxkit::CurveInfo *blowout, const char *saveto);
double dms(const char *pos);
double hms(const char *pos);




#endif


