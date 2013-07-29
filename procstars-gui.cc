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
// Fyi:
// At 8192 px wide, 23 pixels make 1 degree. 1/3 of a pixel is one arc minute, 1/200 of a pixel is 1 arc second.
// At 32768,        91 pixels make 1 degree, 1.5 pixels make one arc minute. 
// A typical human can discerne down to about 1 arcminute.
// the ISS takes up about 30 arc seconds, and is max about -5.9 in magnitude
// The moon takes up about 30 arc minutes, and is about magnitude -12.74
// The sun is about magnitude -26.74, and takes up about 32 arc minutes.
// The dimmest objects observed so far in human visible wavelengths is about magnitude 36.
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
#include <lax/buttondowninfo.h>

#include </usr/include/GraphicsMagick/Magick++.h>

#include <lax/laxoptions.h>
#include <lax/anxapp.h>
#include <lax/curvewindow.h>
//#include <lax/iconselector.h>
#include <lax/button.h>
#include <lax/lineinput.h>
#include <lax/strmanip.h>

#include "catalogs.h"


#define PI (M_PI)
#define TWO_PI (2*M_PI)

//map: scale from [omin,omax] to [nmin,nmax]
#define map(value, omin,omax, nmin,nmax) (nmin+(value-omin)/(omax-omin)*(nmax-nmin))
#define constrain(v, min,max) (v<min ? min : (v>max ? max : v))
#define radians(deg) (deg/180.*M_PI)
#define degrees(rad) (rad*180./M_PI)


#define DBG

using namespace std;
using namespace Magick;
using namespace Laxkit;
using namespace LaxFiles;


//options:


//stuff to make computations easier:
//double pixelwidth=1./width;





//--------------------------------- HaloWindow ---------------------------------
/*! \class HaloWindow
 * Control window for adjusting map of star color index to rgb color.
 */
class HaloWindow : public Laxkit::RowFrame
{
  public:
	RenderContext *context;
	LaxImage *halo;
	CurveWindow *blowout, *ramp;
	ButtonDownInfo buttondown;
	

	HaloWindow(RenderContext *cntxt);
	virtual ~HaloWindow();
	virtual int init();
	virtual void Refresh();
	virtual int Event(const EventData *e,const char *mes);
	virtual void Reset(int which=~0);
	virtual int LBDown(int x,int y,unsigned int state,int count,const LaxMouse *d);
    virtual int LBUp(int x,int y,unsigned int state,const LaxMouse *d);
	virtual int MouseMove(int x,int y,unsigned int state,const LaxMouse *d);

	virtual void UpdateHalo();
	virtual void send();
};

HaloWindow::HaloWindow(RenderContext *cntxt)
  : RowFrame(NULL,"Halo Window","Halo Window",ROWFRAME_ROWS|ANXWIN_DOUBLEBUFFER, 0,310,600,250,0, NULL,0,"update")
{
	context=cntxt;
	halo=create_new_image(150,150);
}

HaloWindow::~HaloWindow()
{
}

int HaloWindow::LBDown(int x,int y,unsigned int state,int count,const LaxMouse *d)
{
	int xx=wholelist.e[0]->x();
	int yy=wholelist.e[0]->y();
	int ww=wholelist.e[0]->w()-2*pad;
	int hh=wholelist.e[0]->h()-2*pad;

	if (x>xx && x<xx+ww && y>yy && y<yy+hh) {
		if (x<xx+ww/2) buttondown.down(d->id, LEFTBUTTON, x,y, -1);
		else buttondown.down(d->id, LEFTBUTTON, x,y, 1);
	}

	return 0;
}

int HaloWindow::LBUp(int x,int y,unsigned int state,const LaxMouse *d)
{
	buttondown.up(d->id, LEFTBUTTON);
	return 0;
}

int HaloWindow::MouseMove(int x,int y,unsigned int state,const LaxMouse *d)
{
	if (!buttondown.any()) return 0;
	int ox,oy;
	buttondown.move(d->id, x,y, &ox,&oy);

	int xx=wholelist.e[0]->x();
	int yy=wholelist.e[0]->y();
	int ww=wholelist.e[0]->w()-2*pad;
	int hh=wholelist.e[0]->h()-2*pad;
	ox-=xx+ww/2;
	oy-=yy+hh/2;
	x-=xx+ww/2;
	y-=yy+hh/2;
	
	double od=sqrt(ox*ox+oy*oy);
	double nd=sqrt(x*x+y*y);

	context->usehalo=(context->usehalo-1)*od/nd + 1;
	if (context->usehalo<1.0001) context->usehalo=1.0001; //otherwise occasional mysterious glitches
	UpdateHalo();
	send();
	needtodraw=1;
	return 0;
}

int HaloWindow::init()
{
	anXWindow *last=NULL;

	CurveWindow *win=NULL;
	
	//AddSpacer(75,25,50,50,  200,100,200,50, -1); //file info
	AddSpacer(200,100,200,50,  200,100,200,50, -1); //preview area


	last=win=ramp=new CurveWindow(NULL,"Ramp","Ramp",0, 5,5,500,500,0, last,object_id,"ramp",
						 			 "Halo Ramp", "radius",0,1,  "a",0,255);
	//radius->AddPoint(.7,255); // *** or whatever halostart is
	ramp->GetInfo()->curvetype=CurveInfo::Autosmooth;
	context->ramp=ramp->GetInfo();
	AddWin(win,1, 200,100,200,50,0,  200,100,200,50,0, -1);


	last=win=blowout=new CurveWindow(NULL,"Blowout","Blowout",0, 5,5,500,500,0, last,object_id,"blowout",
							"Halo Blowout", "opacity",0,1,  "color",0,1);
	context->blowout=blowout->GetInfo();
	blowout->GetInfo()->curvetype=CurveInfo::Autosmooth;
	AddWin(win,1, 200,100,200,50,0,  200,100,200,50,0, -1);

	UpdateHalo();


	last->CloseControlLoop();
	Sync(1);
	return 0;
}

/*! Reset to default colors.
 * which is &1 for red, &2 for green, &4 for blue
 */
void HaloWindow::Reset(int which)
{
	if (which&1) {
		//radius->AddPoint(.7,255); // *** or whatever halostart is
		ramp->Reset();
	}

	if (which&2) {
		blowout->Reset();
	}

	needtodraw=1;
}

void HaloWindow::Refresh()
{
	if (!needtodraw) return;
	if (arrangedstate!=1) Sync(0);

	int pad=10;
	int w=wholelist.e[0]->w()-2*pad;
	int h=wholelist.e[0]->h()-2*pad;
	//double pos, r,g,b;

	if (h<w) w=h;
	foreground_color(0);
	fill_rectangle(this,pad,pad,w,h);
	image_out_skewed(halo, this, pad,pad, w,0, 0,w);

	SwapBuffers();
	needtodraw=0;
}

void HaloWindow::send()
{
    if (win_owner) {
        SimpleMessage *ev=new SimpleMessage;
        app->SendMessage(ev,win_owner,win_sendthis,object_id);
    }
}

void HaloWindow::UpdateHalo()
{
	ramp->GetInfo()->RefreshLookup(256, 0,255);
	blowout->GetInfo()->RefreshLookup(256, 0,255);
	blowout->GetInfo()->LookupDump("blowout",stdout); // ***dbg
	
	unsigned char *data=halo->getImageBuffer();
	int w=halo->w();

	//CreateStockHalo(w, context->usehalo, data, "c", ramp->GetInfo(),blowout->GetInfo(),"halo.png");
	CreateStockHalo(w, context->usehalo, data, "c", ramp->GetInfo(),blowout->GetInfo(),NULL);
	halo->doneWithBuffer(data);

	needtodraw=1;
}

int HaloWindow::Event(const EventData *e,const char *mes)
{
	const SimpleMessage *m=dynamic_cast<const SimpleMessage*>(e);
	if (!m) return anXWindow::Event(e,mes);

	if (!strcmp(mes,"ramp") || !strcmp(mes,"blowout")) {
		UpdateHalo();
		send();
	}

	needtodraw=1;
	return 0;
}




//--------------------------------- IndexWindow ---------------------------------
/*! \class IndexWindow
 * Control window for adjusting map of star color index to rgb color.
 */
class IndexWindow : public Laxkit::RowFrame
{
  public:
	RenderContext *context;
	CurveWindow *rr,*gg,*bb;

	IndexWindow(RenderContext *cntxt);
	virtual ~IndexWindow();
	virtual int init();
	virtual void Refresh();
	virtual int Event(const EventData *e,const char *mes);
	virtual void Reset(int which=~0);
	virtual void UpdateIndex();
	virtual void send();
};

IndexWindow::IndexWindow(RenderContext *cntxt)
  : RowFrame(NULL,"Index Window","Index Window",ROWFRAME_ROWS, 0,0,600,250,0, NULL,0,"update")
{
	context=cntxt;
}

IndexWindow::~IndexWindow()
{
}

int IndexWindow::init()
{
	anXWindow *last=NULL;

	CurveWindow *win=NULL;
	
	AddSpacer(75,25,50,50,  200,150,200,50, -1);


	last=win=rr=new CurveWindow(NULL,"Red","Red",0, 5,5,500,500,0, last,object_id,"red",
						 			 "Red", "B-V",-0.5,2,  "r",0,255);
	context->index_r=rr->GetInfo();
	rr->GetInfo()->curvetype=CurveInfo::Autosmooth;
	rr->AddPoint(.4,255);
	AddWin(win,1, 200,100,200,50,0,  200,150,200,50,0, -1);


	last=win=gg=new CurveWindow(NULL,"Green","Green",0, 5,5,500,500,0, last,object_id,"green",
							"Green", "B-V",-0.5,2,  "g",0,255);
	context->index_g=gg->GetInfo();
	gg->GetInfo()->curvetype=CurveInfo::Autosmooth;
	gg->MovePoint(0, -.5,170);
	gg->MovePoint(1, 2,0);
	gg->AddPoint(.4,255);
	gg->AddPoint(1.7,200);
	AddWin(win,1, 200,100,200,50,0,  200,150,200,50,0, -1);


	last=win=bb=new CurveWindow(NULL,"Blue","Blue",0, 5,5,500,500,0, last,object_id,"blue",
							"Blue", "B-V",-0.5,2, "b",0,255);
	context->index_b=bb->GetInfo();
	bb->GetInfo()->curvetype=CurveInfo::Autosmooth;
	bb->MovePoint(0, -.5,255);
	bb->MovePoint(1, 2,0);
	bb->AddPoint(.5,255);
	bb->AddPoint(1.7,150);
	AddWin(win,1, 200,100,200,50,0,  200,150,200,50,0, -1);

	Reset(~0);

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
		rr->AddPoint(.25,230);
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

	UpdateIndex();

	needtodraw=1;
}

void IndexWindow::Refresh()
{
	if (!needtodraw) return;
	if (arrangedstate!=1) Sync(0);

	int pad=20;
	int w=wholelist.e[0]->w()-2*pad;
	int h=wholelist.e[0]->h()-2*pad;
	double pos, r,g,b;

	StarColor color;
	for (int c=0; c<h; c++) {
		pos=(double)c/(h-1)*(2.5)-.5;

		//------
		indexToRgb(context,pos,0, color);
		r=color.redf();
		g=color.greenf();
		b=color.bluef();
		//------
		//pos=(double)c/(h-1)*255;
		//r=rr->GetInfo()->lookup[(int)pos]/255.;
		//g=gg->GetInfo()->lookup[(int)pos]/255.;
		//b=bb->GetInfo()->lookup[(int)pos]/255.;
		//------
		//r=rr->f(pos)/255;
		//g=gg->f(pos)/255;
		//b=bb->f(pos)/255;

		foreground_color(r,g,b);
		draw_line(this, pad,pad+c, pad+w,pad+c);
	}

	needtodraw=0;
}

void IndexWindow::send()
{
    if (win_owner) {
        SimpleMessage *ev=new SimpleMessage;
        app->SendMessage(ev,win_owner,win_sendthis,object_id);
    }
}

void IndexWindow::UpdateIndex()
{
	rr->GetInfo()->RefreshLookup(256, 0,255);
	gg->GetInfo()->RefreshLookup(256, 0,255);
	bb->GetInfo()->RefreshLookup(256, 0,255);
}

int IndexWindow::Event(const EventData *e,const char *mes)
{
	const SimpleMessage *m=dynamic_cast<const SimpleMessage*>(e);
	if (!m) return anXWindow::Event(e,mes);

	if (!strcmp(mes,"red") || !strcmp(mes,"green") || !strcmp(mes,"blue")) {
		send();
	}

	needtodraw=1;
	return 0;
}


//--------------------------------- SizeWindow ---------------------------------
/*! \class SizeWindow
 * Control window for adjusting sizes of stars.
 */
class SizeWindow : public Laxkit::RowFrame
{
  public:
	RenderContext *context;

	CurveWindow *bigscale,*pointopacity;
	SizeWindow(RenderContext *rr);
	virtual ~SizeWindow();
	virtual int init();
	virtual void Refresh();
	virtual int Event(const EventData *e,const char *mes);
	virtual void Reset(int which=~0);
	virtual void UpdateSize();
	virtual void send(const char *mes=NULL);
};

SizeWindow::SizeWindow(RenderContext *rr)
  : RowFrame(NULL,"Size Window","Size Window",ROWFRAME_ROWS, 0,0,600,250,0, NULL,0,"update")
{
	context=rr;
	bigscale=NULL;
	pointopacity=NULL;
}

SizeWindow::~SizeWindow()
{
}

int SizeWindow::init()
{
	anXWindow *last=NULL;

	CurveWindow *win=NULL;
	
	AddSpacer(75,25,50,50,  200,150,200,50, -1); //*** for range selector


	last=win=bigscale=new CurveWindow(NULL,"Large","Large",0, 5,5,500,500,0, last,object_id,"large",
						 			 "Big star scale", "magnitude",context->bigthreshhold,-1, "px",1,100);
	bigscale->editable=CurveWindow::YMax|CurveWindow::YUnits;
	bigscale->GetInfo()->MovePoint(1, context->maximum_magnitude,context->maxstarsize);
	context->bigscale=bigscale->GetInfo();
	win->GetInfo()->curvetype=CurveInfo::Autosmooth;
	AddWin(win,1, 200,100,200,50,0,  200,150,200,50,0, -1);


	last=win=pointopacity=new CurveWindow(NULL,"Points","Points",0, 5,5,500,500,0, last,object_id,"points",
							"Point Opacity", "mag",context->minimum_magnitude,context->bigthreshhold,  "a",0,1);
	context->pointopacity=pointopacity->GetInfo();
	win->GetInfo()->curvetype=CurveInfo::Autosmooth;
	AddWin(win,1, 200,100,200,50,0,  200,150,200,50,0, -1);


//	last=win=pointamp=new CurveWindow(NULL,"Point Amp","Point Amp",0, 5,5,500,500,0, last,object_id,"pointamp",
//							"Point Amp", "mag",context->minimum_magnitude,context->bigthreshhold, "a",0,1);
//	context->pointamp=pointamp->GetInfo();
//	win->GetInfo()->curvetype=CurveInfo::Autosmooth;
//	win->MovePoint(1, 0,0);
//	AddWin(win,1, 200,100,200,50,0,  200,150,200,50,0, -1);

	UpdateSize();

	last->CloseControlLoop();
	Sync(1);
	return 0;
}

/*! Reset to default colors.
 * which is &1 for red, &2 for green, &4 for blue
 */
void SizeWindow::Reset(int which)
{
	if (which&1) {
	}
	needtodraw=1;
}

void SizeWindow::Refresh()
{
	if (!needtodraw) return;
	if (arrangedstate!=1) Sync(0);

	needtodraw=0;
}

void SizeWindow::send(const char *mes)
{
    if (win_owner) {
        SimpleMessage *ev=new SimpleMessage;
        app->SendMessage(ev,win_owner,mes?mes:win_sendthis,object_id);
    }
}

void SizeWindow::UpdateSize()
{
	bigscale    ->GetInfo()->RefreshLookup(256, 1,bigscale->GetInfo()->ymax);
	pointopacity->GetInfo()->RefreshLookup(256, 0,255);
}

int SizeWindow::Event(const EventData *e,const char *mes)
{
	const SimpleMessage *m=dynamic_cast<const SimpleMessage*>(e);
	if (!m) return anXWindow::Event(e,mes);

	if (!strcmp(mes,"large") || !strcmp(mes,"points")) {
		UpdateSize();
		send();
	}

	if (!strcmp(mes,"large")) {
		send("large");
	}


	needtodraw=1;
	return 0;
}




//--------------------------------- MainWindow ---------------------------------
/*! \class MainWindow
 */
class MainWindow : public Laxkit::RowFrame
{
  public:
	int firsttime;

	int numstars;
	RandomCatalog previewcatalog;
	LaxImage *preview;
	int needtoupdate;

	HaloWindow *halowindow;

	RenderContext *context;
	MainWindow(RenderContext *rr);
	virtual ~MainWindow() {}
	virtual int init();
	virtual void Refresh();
	virtual int Event(const EventData *e,const char *mes);
	//virtual void Reset(int which=~0);
	virtual void UpdatePreview();

    virtual int LBUp(int x,int y,unsigned int state,const LaxMouse *d);
};

MainWindow::MainWindow(RenderContext *rr)
  : RowFrame(NULL,"Main Window","Main Window",ANXWIN_ESCAPABLE|ROWFRAME_ROWS, 500,0,600,900,0, NULL,0,NULL),
	previewcatalog("preview",1000)
{
	numstars=1000;
	context=rr;
	preview=NULL;
	firsttime=2;
	needtoupdate=1;
}

int MainWindow::LBUp(int x,int y,unsigned int state,const LaxMouse *d)
{
	if (x<win_w/3) numstars*=.6;
	else if (x>win_w*2/3) numstars*=1.3;
	if (numstars<100) numstars=100;

	cerr <<"Regenerate "<<numstars<<"stars..."<<endl;


	previewcatalog.Repopulate(numstars,0);
	UpdatePreview();
	needtodraw=1;
	return 0;
}

int MainWindow::init()
{
	anXWindow *last=NULL;

	AddSpacer(1000,2005,1000,50,  400,200,200,50, -1); //area for preview image
	AddNull();


	//  Save To: ____ [...]    <<RENDER NOW>>    Save Settings To: ______[...] [load]
	LineInput *file;
	last=file=new LineInput(this, "File","File",0, 0,0,0,0,0, NULL,object_id,"file","Render To:","stars#.tif");
	AddWin(file,1, file->win_w,0,300,50,0, file->win_h,0,10,50,0, -1);

	Button *button;
	last=button=new Button(this, "Render","Render",0, 0,0,0,0,0, last,object_id,"render",0, " Render now! ");
	button->gap=button->win_w*.25;
	AddWin(button,1, button->win_w,0,300,50,0, button->win_h,0,10,50,0, -1);

	last=file=new LineInput(this, "Save Settings","Save Settings",0, 0,0,0,0,0, last,object_id,"settings","Project:","");
	AddWin(file,1, file->win_w,0,300,50,0, file->win_h,0,10,50,0, -1);

	AddNull();


	 //-----editing windows:

	SizeWindow *swin=new SizeWindow(context);
	swin->SetOwner(this);
	AddWin(swin,1, 800,400,200,50,0,  200,150,200,50,0, -1);

	IndexWindow *iwin=new IndexWindow(context);
	iwin->SetOwner(this);
	AddWin(iwin,1, 800,400,200,50,0,  200,150,200,50,0, -1);

	halowindow=new HaloWindow(context);
	halowindow->SetOwner(this);
	AddWin(halowindow,1, 800,400,200,50,0,  200,150,200,50,0, -1);

	//CatalogWindow *catwin=new CatalogWindow();
	//AddWin(catwin,1, 800,400,200,50,0,  200,150,200,50,0, -1);



	last->CloseControlLoop();
	Sync(1);
	return 0;
}

int MainWindow::Event(const EventData *e,const char *mes)
{
	const SimpleMessage *m=dynamic_cast<const SimpleMessage*>(e);
	if (!m) return anXWindow::Event(e,mes);

	if (!strcmp(mes,"render")) {
		cout <<"Render..."<<endl;
		return 0;
	}

	if (!strcmp(mes,"large")) {
		 //max star size changed, need to update halo sample
		//halowindow->***;
	}

	if (!strcmp(mes,"update")) {
		needtoupdate=1;
		needtodraw=1;
		return 0;
	}

	return anXWindow::Event(e,mes);
}

void MainWindow::UpdatePreview()
{
	cerr <<"Updating preview..."<<endl;

	if (context->halo && context->halowidth<context->maxstarsize*context->usehalo) {
		 //star size was changed, we need to make a different size halo
		cerr << "Resizing halo sample size.."<<endl;
		delete[] context->halo;
		context->halo=NULL;
	}
	if (!context->halo) {
		halowindow->UpdateHalo();
		if (!context->halo) {
			context->halowidth=2*context->maxstarsize*context->usehalo;
			context->halo=new unsigned char[context->halowidth*context->halowidth];
		}
	}
	CreateStockHalo(context->halowidth, context->usehalo, context->halo, "g",
					context->ramp,context->blowout,"halo.png");

	int ww=wholelist.e[0]->w();
	int hh=wholelist.e[0]->h();

	if (preview && (ww!=preview->w() || hh!=preview->h())) {
		preview->dec_count();
		preview=NULL;
	}
	if (!preview) {
		preview=create_new_image(ww,hh);
	}

	unsigned char *data=preview->getImageBuffer();
	ww=preview->w();
	hh=preview->h();
	memset(data,0,ww*hh*4);
	for (int x=0; x<ww; x++) {
	  for (int y=0; y<hh; y++) {
		data[(y*ww+x)*4+3]=255;
	  }
	}

	previewcatalog.Render(context, data,ww,hh);
	preview->doneWithBuffer(data);
}

void MainWindow::Refresh()
{
	if (!needtodraw) return;
	if (arrangedstate!=1) Sync(0);

	if (firsttime>0) {
		firsttime--;
		return;
	}


	needtodraw=0;

	int xx=wholelist.e[0]->x();
	int yy=wholelist.e[0]->y();
	int ww=wholelist.e[0]->w();
	int hh=wholelist.e[0]->h();

	if (!preview || needtoupdate) {
		UpdatePreview();
		needtoupdate=0;
	}

	image_out_skewed(preview, this, xx,yy, ww,0, 0,hh);
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
	double alphaamp=0;


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
				alphaamp=strtod(o->arg(),NULL);
				cerr <<" *** need to implement putting on base alphaamp to point amp curve:"<<alphaamp<<endl;
              } break;

		}
	}


	if (ang>0) {
		rr.maxstarsize=rr.width*ang/60./360.;
	}
	if (rr.usehalo>1) {
		rr.halowidth=rr.maxstarsize*rr.usehalo;
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
	makestr(app.app_profile,"Dark");
	app.init(argc,argv);

	app.addwindow(new MainWindow(&rr));
	

	 //...and off we go!
	app.run();


	cout <<"--- Shutting down... ----\n";
	app.close();
	
	cout <<"------ Bye! -------\n";
	return 0;

} //main()










