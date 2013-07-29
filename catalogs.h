


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
	double red() { return r; }
	void red(double rr) { r=rr; }
	double green() { return g; }
	void green(double gg) { g=gg; }
	double blue() { return b; }
	void blue(double bb) { b=bb; }
	double alpha() { return a; }
	void alpha(double aa) { a=aa; }
};

//------------------------------- RenderContext -----------------------------------
class Catalog;

class RenderContext
{
  public:
	char *projectfile; //save settings here
	char *filename; //final stars file
	long width;
	long height;
	unsigned char *data;

	int galactic;
	int transparent;
	double magnitude;
	double maxmagnitude;
	double maxstarsize; //pixels of biggest star
	double bigthreshhold; //magnitude < this get big treatment
	double alphaamp;

	double usehalo;
	unsigned char *halo;
	int halowidth; //width in pixels of halo image

	Laxkit::CurveInfo *ramp;
	Laxkit::CurveInfo *blowout;
	Laxkit::CurveInfo *index_r;
	Laxkit::CurveInfo *index_g;
	Laxkit::CurveInfo *index_b;

	double min_asc;
	double max_asc;
	double min_dec;
	double max_dec;


	double minimum_magnitude;
	double maximum_magnitude;

	Laxkit::PtrStack<Catalog> catalogs;
	Catalog *catalog; //current

	RenderContext();
};



//------------------------------- Catalog -----------------------------------
enum CatalogTypes
{
	Tycho2,
	PrincipalGalaxy,
	RandomMemory,
	CustomStars,
	CustomGalaxies
};

class Catalog
{
  protected:
  public:
	char *name;
	char *filename;
	CatalogTypes type;
	
	int *magnitude_distribution;
	int nummags;
	double minimum_magnitude;
	double maximum_magnitude;

	int min_mag_cutoff;
	int max_mag_cutoff;

	Catalog(const char *nname, const char *nfile, CatalogTypes ntype);
	virtual ~Catalog();

	virtual int Render(RenderContext *context);
	//virtual int GetStats() = 0;

	virtual int OpenCatalog();
	//virtual int GetLine(int &object_type, double &asc, double &dec, double &bmag, double &vmag) = 0;
	virtual int CloseCatalog();
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

flatvector Eq2Gal(float ra, float dec);
void indexToRgb(RenderContext *rr, float index, float vmag, StarColor &color);
double indexToRed(float index);
double indexToGreen(float index);
double indexToBlue(float index);
void drawStar(RenderContext *rr, float ra, float dec, float vmag, float bmag);
void drawStarSimple(RenderContext *rr, double ra, double dec, double index, double vmag);
void CreateStockHalo(int w,double halosize, unsigned char *halo, const char *format, Laxkit::CurveInfo *ramp, Laxkit::CurveInfo *blowout);
double dms(const char *pos);
double hms(const char *pos);

int Render(RenderContext *context);
int Process_PGC(RenderContext *context);
int Process_Tycho(RenderContext *context);



#endif


