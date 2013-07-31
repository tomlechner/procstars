


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

	RenderContext();
	virtual ~RenderContext();
	virtual void Reset(int which);

	virtual void dump_out(FILE *f,int indent,int what,Laxkit::anObject *context);
    virtual LaxFiles::Attribute *dump_out_atts(LaxFiles::Attribute *att,int what,Laxkit::anObject *context);
    virtual void dump_in_atts(LaxFiles::Attribute *att,int flag,Laxkit::anObject *context);
};



//------------------------------- Catalog -----------------------------------
enum CatalogTypes
{
	Unknown,
	Tycho2,
	PrincipalGalaxy,
	RandomMemory,
	CustomStars,
	CustomGalaxy
};

class Catalog
{
  protected:
  public:
	char *name;
	char *filename;
	CatalogTypes type;
	
	int visible;

	 //stats
	double minimum_magnitude;
	double maximum_magnitude;
	int numstars;

	 //subset
	double min_mag_cutoff;
	double max_mag_cutoff;

	Catalog(const char *nname, const char *nfile, CatalogTypes ntype);
	virtual ~Catalog();

	virtual const char *TypeName();
	virtual int Render(RenderContext *context);
	//virtual int GetStats() = 0;

	virtual int OpenCatalog();
	//virtual int GetLine(int &object_type, double &asc, double &dec, double &bmag, double &vmag) = 0;
	virtual int CloseCatalog();

	virtual void dump_out(FILE *f,int indent,int what,Laxkit::anObject *context);
    virtual LaxFiles::Attribute *dump_out_atts(LaxFiles::Attribute *att,int what,Laxkit::anObject *context);
    virtual void dump_in_atts(LaxFiles::Attribute *att,int flag,Laxkit::anObject *context);
};

class RandomCatalog : public Catalog
{
  public:
	int numpoints;
	double *color_index;
	double *color_mag;
	double *asc;
	double *dec;

	RandomCatalog(const char *nname, int num);
	virtual int Render(RenderContext *context);
	virtual int Render(RenderContext *context, unsigned char *data,int ww,int hh);

	virtual int Repopulate(int num, int spherical);
};

//------------------------------- Rendering Misc -----------------------------------

flatvector Eq2Gal(double ra, double dec);
void indexToRgb(RenderContext *rr, double index, double vmag, StarColor &color);
void drawStar(RenderContext *rr, double ra, double dec, double vmag, double bmag);
void drawStarSimple(RenderContext *context, double ra, double dec, double index, double vmag);
void CreateStockHalo(int w,double halosize, unsigned char *halo, const char *format,
					 Laxkit::CurveInfo *ramp, Laxkit::CurveInfo *blowout, const char *saveto);
double dms(const char *pos);
double hms(const char *pos);

int Render(RenderContext *context);
int Process_Galaxy(RenderContext *rr);
int Process_PGC(RenderContext *context);
int Process_Tycho(RenderContext *context);



#endif


